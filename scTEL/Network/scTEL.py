from numpy import inf, zeros, zeros_like as np_zeros_like, empty
from pandas import concat
from anndata import AnnData
from torch import cat, no_grad, randn, zeros_like, zeros as torch_zeros, argmax, save
from torch.nn import Module, Linear, Sequential, Softmax, LSTMCell
from torch.optim import Adam
from torch.optim.lr_scheduler import StepLR

from .TE_Layer import Input_TE, FF_TE, LambdaLayer, Dual_Forward, Protein_Linear


class scTEL_Model(Module):
    def __init__(self, p_mod1, p_mod2, h_size, h, loss1, loss2, quantiles, categories, drop_rate=0.2):
        # p_mod1: HVGs p_mod2:ref_proteins
        super(scTEL_Model, self).__init__()

        h_size, drop_rate = h_size, drop_rate

        self.LSTMCell1 = LSTMCell(h_size, h_size)
        self.LSTMCell2 = LSTMCell(h_size, h_size)
        self.LSTMCell3 = LSTMCell(h_size, h_size)

        self.input_TE_1 = Input_TE(p_mod1, h_size, h, drop_rate)
        self.FF_TE_2 = FF_TE(h_size, h, drop_rate)
        self.FF_TE_3 = FF_TE(h_size, h, drop_rate)

        MSE_output = Protein_Linear(h_size, p_mod2, dropout_rate=drop_rate)

        if len(quantiles) > 0:
            quantile_layer = []
            quantile_layer.append(Linear(h_size, p_mod2 * len(quantiles)))  # 512→224×4
            quantile_layer.append(LambdaLayer(lambda x: x.view(-1, p_mod2, len(quantiles))))  # 224×4→(-1,224,4)
            quantile_layer = Sequential(*quantile_layer)

            self.mod2_out = Dual_Forward(MSE_output, quantile_layer)

        else:
            self.mod2_out = MSE_output

        if categories is not None:
            self.celltype_out = Sequential(Linear(h_size, len(categories)), Softmax(1))
            self.forward = self.forward_transfer

            self.categories_arr = empty((len(categories),), dtype='object')
            for cat in categories:
                self.categories_arr[categories[cat]] = cat

        else:
            self.forward = self.forward_simple
            self.categories_arr = None

        self.quantiles = quantiles
        self.loss1, self.loss2 = loss1, loss2

    def forward_transfer(self, x):
        x = self.input_TE_1(x)
        h, c = self.LSTMCell1(x, (zeros_like(x), zeros_like(x)))

        x = self.FF_TE_2(x)
        h, c = self.LSTMCell2(x, (h, c))

        x = self.FF_TE_3(x)
        h, c = self.LSTMCell3(x, (h, c))



        return {'celltypes': self.celltype_out(c), 'modality 2': self.mod2_out(x), 'embedding': h}

    def forward_simple(self, x):
        x = self.input_TE_1(x)
        h, c = self.LSTMCell1(x, (zeros_like(x), zeros_like(x)))

        x = self.FF_TE_2(x)
        h, c = self.LSTMCell2(x, (h, c))

        x = self.FF_TE_3(x)
        h, c = self.LSTMCell3(x, (h, c))
        

        return {'celltypes': None, 'modality 2': self.mod2_out(x), 'embedding': h}

    def train_backprop(self, train_loader, val_loader, path,
                       n_epoch=1000, ES_max=30, decay_max=6, decay_step=0.1, lr=10 ** (-3)):
        optimizer = Adam(self.parameters(), lr=lr)
        scheduler = StepLR(optimizer, step_size=1, gamma=decay_step)

        patience = 0
        bestloss = inf

        if self.categories_arr is None:
            get_correct = lambda x: 0
        else:
            get_correct = lambda outputs: (argmax(outputs['celltypes'], axis=1) == celltypes).sum()

        for epoch in range(n_epoch):
            with no_grad():
                running_loss, rtype_acc, protein_loss, type_loss = 0., 0., 0., 0.
                self.eval()
                for batch, inputs in enumerate(val_loader):
                    mod1, mod2, protein_bools, celltypes = inputs
                    outputs = self(mod1)

                    n_correct = get_correct(outputs)
                    mod1_loss = self.loss1(outputs['celltypes'], celltypes)
                    mod2_loss = self.loss2(outputs['modality 2'], mod2, protein_bools)

                    rtype_acc += n_correct
                    type_loss  += mod1_loss.item() * len(mod1)
                    protein_loss += mod2_loss.item() * len(mod2)
                    running_loss += (mod1_loss + mod2_loss).item() * len(mod2)

                if self.categories_arr is None:
                    print(
                        f"Epoch {epoch} val total loss = {running_loss / len(val_loader):.3f}")
                else:
                    print(f"Epoch {epoch} val total loss = {running_loss / len(val_loader):.3f},Loss type = {type_loss / len(val_loader):.3f},"
                          f"Loss protein = {protein_loss / len(val_loader):.3f}, validation accuracy = {rtype_acc / len(val_loader):.3f}")

                patience += 1

                if bestloss > running_loss:
                    bestloss, patience = running_loss, 0
                    save(self.state_dict(), path)

                if (patience + 1) % decay_max == 0:
                    scheduler.step()
                    print(f"Decaying loss to {optimizer.param_groups[0]['lr']}")

                if (patience + 1) > ES_max:
                    break

            self.train()
            train_loss, train_type_loss, train_protein_loss,rtype_acc = 0., 0., 0.,0.
            for batch, inputs in enumerate(train_loader):
                optimizer.zero_grad()

                mod1, mod2, protein_bools, celltypes = inputs
                outputs = self(mod1)
                n_correct = get_correct(outputs)
                rtype_acc += n_correct
                mod1_loss = self.loss1(outputs['celltypes'], celltypes)
                mod2_loss = self.loss2(outputs['modality 2'], mod2, protein_bools)

                train_type_loss += mod1_loss.item() * len(mod1)
                train_protein_loss += mod2_loss.item() * len(mod2)
                loss = mod1_loss + mod2_loss
                train_loss += loss.item() * len(mod2)
                loss.backward()

                optimizer.step()
            if self.categories_arr is None:
                print(f"Epoch {epoch} train total loss = {train_protein_loss / len(train_loader):.3f}")
            else:
                print(f"Epoch {epoch} train total loss = {train_loss / len(train_loader):.3f}, Loss type = {train_type_loss / len(train_loader):.3f},"
                      f" Loss protein = {train_protein_loss / len(train_loader):.3f}, train accuracy = {rtype_acc / len(train_loader):.3f}")

    def impute(self, impute_loader, requested_quantiles, proteins):
        imputed_test = proteins.copy()

        for quantile in requested_quantiles:
            imputed_test.layers['q' + str(round(100 * quantile))] = np_zeros_like(imputed_test.X)

        self.eval()
        start = 0

        for mod1, bools, celltypes in impute_loader:
            end = start + mod1.shape[0]

            with no_grad():
                outputs = self(mod1)

            if len(self.quantiles) > 0:
                mod2_impute, mod2_quantile = outputs['modality 2']
            else:
                mod2_impute = outputs['modality 2']

            imputed_test.X[start:end] = self.fill_predicted(imputed_test.X[start:end], mod2_impute, bools)

            for quantile in requested_quantiles:
                index = [i for i, q in enumerate(self.quantiles) if quantile == q][0]
                q_name = 'q' + str(round(100 * quantile))
                imputed_test.layers[q_name][start:end] = mod2_quantile[:, :, index].cpu().numpy()

            start = end

        return imputed_test

    def embed(self, impute_loader, test_loader, cells_train, cells_test):
        if cells_test is not None:
            embedding = AnnData(zeros(shape=(len(cells_train) + len(cells_test), 512)))
            embedding.obs = concat((cells_train, cells_test), join='inner')
        else:
            embedding = AnnData(zeros(shape=(len(cells_train), 512)))
            embedding.obs = cells_train

        self.eval()
        start = 0

        for mod1, bools, celltypes in impute_loader:
            end = start + mod1.shape[0]
            outputs = self(mod1)

            embedding[start:end] = outputs['embedding'].detach().cpu().numpy()

            start = end

        if cells_test is not None:
            for mod1 in test_loader:
                end = start + mod1.shape[0]
                outputs = self(mod1)

                embedding[start:end] = outputs['embedding'].detach().cpu().numpy()

                start = end

        return embedding

    def fill_predicted(self, array, predicted, bools):
        bools = bools.cpu().numpy()
        return (1. - bools) * predicted.cpu().numpy() + array

    def predict(self, test_loader, requested_quantiles, proteins, cells):
        imputed_test = AnnData(zeros(shape=(len(cells), len(proteins.var))))
        imputed_test.obs = cells
        imputed_test.var.index = proteins.var.index

        if self.categories_arr is not None:
            celltypes = ['None'] * len(cells)

        for quantile in requested_quantiles:
            imputed_test.layers['q' + str(round(100 * quantile))] = np_zeros_like(imputed_test.X)

        self.eval()
        start = 0

        for mod1 in test_loader:
            end = start + mod1.shape[0]

            with no_grad():
                outputs = self(mod1)

            if self.categories_arr is not None:
                predicted_types = argmax(outputs['celltypes'], axis=1).cpu().numpy()
                celltypes[start:end] = self.categories_arr[predicted_types].tolist()

            if len(self.quantiles) > 0:
                mod2_impute, mod2_quantile = outputs['modality 2']
            else:
                mod2_impute = outputs['modality 2']

            imputed_test.X[start:end] = mod2_impute.cpu().numpy()

            for quantile in requested_quantiles:
                index = [i for i, q in enumerate(self.quantiles) if quantile == q][0]
                q_name = 'q' + str(round(100 * quantile))
                imputed_test.layers[q_name][start:end] = mod2_quantile[:, :, index].cpu().numpy()

            start = end

        if self.categories_arr is not None:
            imputed_test.obs['transfered cell labels'] = celltypes

        return imputed_test
