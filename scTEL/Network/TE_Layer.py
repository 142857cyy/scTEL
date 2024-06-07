from torch import cat
from torch.nn import Module, Linear, BatchNorm1d, PReLU, Dropout, TransformerEncoderLayer


class Input_Embedding(Module):
    def __init__(self, in_units, out_units, dropout=0.15): 
        super(Input_Embedding, self).__init__()

        self.embedding = Linear(in_units, out_units)
        self.bnorm_in = BatchNorm1d(out_units)
        self.Linear_out = Linear(1, out_units)
        self.bnorm_out = BatchNorm1d(out_units)
        self.act = PReLU()
        self.dropout_out = Dropout(dropout)

    def forward(self, x):
        x = self.embedding(x)  
        x = self.bnorm_in(x)
        x = self.act(x)
        x = x.unsqueeze(-1)  
        x = self.Linear_out(x)  
        x = self.bnorm_out(x)
        x = self.act(x)
        x = self.dropout_out(x)

        return x


class Input_TE(Module):
    def __init__(self, in_units, out_units, h, dropout_rate):  
        super(Input_TE, self).__init__()

        self.embedding = Input_Embedding(in_units, out_units, dropout=dropout_rate)
        self.bnorm_out = BatchNorm1d(out_units)
        self.act = PReLU(num_parameters=out_units)
        self.dropout_out = Dropout(dropout_rate)
        self.Linear_in = Linear(1, out_units)
        self.TE = TransformerEncoderLayer(d_model=out_units, nhead=h, dim_feedforward=out_units, dropout=0.2,
                                          batch_first=True)
        self.Linear_out = Linear(out_units, 1)

    def forward(self, x):
        x = self.embedding(x)
        x = self.TE(x)  
        x = self.Linear_out(x)  
        x = self.bnorm_out(x)
        x = self.act(x)
        x = x.squeeze(-1) 

        return x


class FF_Embedding(Module):
    def __init__(self, d):
        super(FF_Embedding, self).__init__()

        self.FC = Linear(1, d)
        self.BN = BatchNorm1d(d)
        self.act = PReLU()

    def forward(self, x):
        x = x.unsqueeze(-1)
        x = self.FC(x)
        x = self.act(x)

        return x


class FF_TE(Module):
    def __init__(self, hidden_units, h, dropout_rate):
        super(FF_TE, self).__init__()

        self.FF_Em = FF_Embedding(d=hidden_units)
        self.TE = TransformerEncoderLayer(d_model=hidden_units, nhead=h, dim_feedforward=hidden_units, dropout=0.2,
                                          batch_first=True)
        self.FC = Linear(hidden_units, 1)
        self.BN = BatchNorm1d(hidden_units)
        self.act = PReLU(num_parameters=hidden_units)
        self.Dropout = Dropout(dropout_rate)

    def forward(self, x):
        x = self.FF_Em(x)
        x = self.TE(x)  
        x = self.FC(x)  
        x = self.BN(x)
        x = self.act(x)
        x = self.Dropout(x)
        x = x.squeeze(-1) 

        return x


class Protein_Linear(Module):
    def __init__(self, d, p_mod2, dropout_rate=0.15):
        super(Protein_Linear, self).__init__()

        self.FC1 = Linear(d, d)
        self.BN1 = BatchNorm1d(d)
        self.FC2 = Linear(d, 256)
        self.BN2 = BatchNorm1d(256)
        self.Dropout = Dropout(dropout_rate)
        self.FC3 = Linear(256, p_mod2)
        self.act = PReLU()

    def forward(self, x):
        x = self.FC1(x)
        x = self.BN1(x)
        x = self.act(x)
        x = self.Dropout(x)
        x = self.FC2(x)
        x = self.act(x)
        x = self.FC3(x)
        x = self.act(x)

        return x


class LambdaLayer(Module):
    def __init__(self, lambd):
        super(LambdaLayer, self).__init__()
        self.lambd = lambd

    def forward(self, x):
        return self.lambd(x)


class Dual_Forward(Module):
    def __init__(self, MSE_layer, Quantile_layer):
        super(Dual_Forward, self).__init__()
        self.MSE_layer = MSE_layer
        self.Quantile_layer = Quantile_layer

    def forward(self, decoded):
        return self.MSE_layer(decoded), self.Quantile_layer(decoded)
