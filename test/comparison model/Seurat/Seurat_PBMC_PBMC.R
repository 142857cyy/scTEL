library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(SeuratDisk)


reference = LoadH5Seurat("../Data/pbmc_multimodal.h5seurat")  # (161159,20729)
reference = subset(reference, celltype.l2 != 'Doublet')  
donor = reference@meta.data['donor']

train_donors = c('P1','P3','P4','P7')
bools = sapply(donor, function(x) x %in% train_donors)[,1]

gene_matrix = reference@assays$SCT@counts[,bools]
protein_matrix = reference@assays$ADT@counts[,bools]

# 筛选基因
cell_indices = colSums(gene_matrix >= 1) >= 200
gene_matrix = gene_matrix[, cell_indices]  # (20729,75483)
protein_matrix = protein_matrix[, cell_indices] # (224,75483)

# 提取训练集的细胞meta信息
meta_data = reference@meta.data[bools,][cell_indices,] # (75483,12)
# 创建训练集Seurat对象
train_assay = CreateSeuratObject(
  gene_matrix,
  assay = 'RNA',
  meta.data = meta_data
)
protein =CreateSeuratObject(
  protein_matrix,
  assay = 'ADT',
  meta.data = meta_data
)

# 筛选基因/蛋白质
keep_features = rownames(train_assay@assays$RNA$counts)[rowSums(train_assay@assays$RNA$counts >= 1) >= 30]
train_assay = subset(train_assay, features = keep_features) # (17637,75483)
keep_features = rownames(protein@assays$ADT$counts)[rowSums(protein@assays$ADT$counts) >= 1]
protein = subset(protein, features = keep_features)  # (224,75483)
train_assay@assays$ADT = protein@assays$ADT
train_list <- SplitObject(train_assay, split.by = 'donor')

rm(keep_features)
rm(meta_data)
rm(cell_indices)
rm(gene_matrix)
rm(protein_matrix)
rm(train_assay)
rm(protein)

for(i in 1:length(train_list)){
  DefaultAssay(train_list[[i]]) <- 'RNA'
  train_list[[i]] <- SCTransform(train_list[[i]], verbose =TRUE)
  
  DefaultAssay(train_list[[i]]) <- 'ADT'
  VariableFeatures(train_list[[i]]) <- rownames(train_list[[i]][['ADT']])
  train_list[[i]] <- NormalizeData(train_list[[i]], normalization.method = 'CLR', margin = 2) %>% ScaleData()
}

# Integrate RNA
for (i in 1:length(train_list)) {
  DefaultAssay(train_list[[i]]) <- 'SCT'
}

train_int_features <- SelectIntegrationFeatures(
  object.list = train_list, nfeatures = 3000)
options(future.globals.maxSize= 30*1024^3)
train_list <- PrepSCTIntegration(
  object.list = train_list, anchor.features = train_int_features, 
  verbose = FALSE)
train_list <- lapply(X = train_list, FUN = function(x) {
  x <- RunPCA(x, features = train_int_features, verbose = FALSE)
})

train_anchors <- FindIntegrationAnchors(
  object.list = train_list, reduction = "rpca",
  anchor.features = train_int_features,
  scale=FALSE, normalization.method = 'SCT',
  dims = 1:30)
train_rna_integrated <- IntegrateData(anchorset = train_anchors, dims = 1:30)

# Integrate Protein
for (i in 1:length(train_list)) {
  DefaultAssay(train_list[[i]]) <- 'ADT'
}
train_list <- lapply(X = train_list, FUN = function(x) {
  x <- RunPCA(x, reduction.name = 'apca', verbose = FALSE)
})
for (i in 1:length(train_list)) {
  train_list[[i]]@reductions$rpca <- train_list[[i]]@reductions$pca
  train_list[[i]]@reductions$pca <- train_list[[i]]@reductions$apca
}
prt_int_features <- rownames(train_list[[1]]@assays$ADT)
prt_anchors <- FindIntegrationAnchors(
  object.list = train_list, reduction = "rpca",
  anchor.features = prt_int_features,
  scale=FALSE, normalization.method = 'LogNormalize',
  dims = 1:30)

rm(train_int_features)
rm(prt_int_features)
rm(train_list)
train_protein_integrated <- IntegrateData(anchorset = prt_anchors, dims = 1:30)

# Combine Integrations
train_assay <- CreateSeuratObject(
  counts = train_rna_integrated@assays$integrated,
  assay = 'RNA',
  meta.data = train_rna_integrated@meta.data,
  min.cells=0,
  min.features=0)
train_assay@assays$ADT <- train_protein_integrated@assays$integrated
DefaultAssay(train_assay) = "RNA"

# Create Reference Data
train_assay = ScaleData(train_assay) %>% RunPCA()
DefaultAssay(train_assay) <- 'ADT'
train_assay <- ScaleData(train_assay) %>% RunPCA(reduction.name = 'apca')
train_assay <- FindMultiModalNeighbors(
  train_assay, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)
train_assay <- RunUMAP(train_assay, nn.name = "weighted.nn", return.model=T,
                       reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

train_assay <- FindClusters(train_assay, graph.name = "wsnn",
                            algorithm = 3, resolution = 0.2, verbose = FALSE)

options(repr.plot.width=7, repr.plot.height=7)

p1 <- DimPlot(train_assay, reduction = 'wnn.umap', group.by = "orig.ident")
p1

# Run supervised PCA on weighted NN graph
DefaultAssay(train_assay) <- 'RNA'
train_assay <- ScaleData(train_assay, assay = 'RNA')
train_assay <- RunSPCA(train_assay, assay='RNA', graph='wsnn')

train_assay <- FindNeighbors(
  object = train_assay,
  reduction = "spca",
  dims = 1:50,
  graph.name = "spca.annoy.neighbors", 
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)
rm(train_rna_integrated)
rm(train_protein_integrated)

# Imoute Query Dataset
gene_matrix = reference@assays$SCT@counts[, !bools]
protein_matrix = reference@assays$ADT@counts[, !bools]
cell_indices = colSums(gene_matrix >= 1) >= 200
gene_matrix = gene_matrix[, cell_indices]
protein_matrix = protein_matrix[, cell_indices]
print(dim(gene_matrix))
meta_data = reference@meta.data[!bools,][cell_indices,]

test_assay = CreateSeuratObject(
  gene_matrix,
  assay = "RNA",
  meta.data = meta_data)
rm(gene_matrix)
protein = CreateSeuratObject(
  protein_matrix,
  assay = "ADT",
  meta.data = meta_data)
rm(protein_matrix)

keep_features = rownames(test_assay@assays$RNA$counts)[rowSums(test_assay@assays$RNA$counts >= 1) >= 30]
test_assay = subset(test_assay, features = keep_features)
keep_features = rownames(protein@assays$ADT$counts)[rowSums(protein@assays$ADT$counts) >= 1]
protein = subset(protein, features = keep_features)
test_assay@assays$ADT = protein@assays$ADT #comment this
test_list <- SplitObject(test_assay, split.by = "donor")
for (i in 1:length(test_list)) {
  DefaultAssay(test_list[[i]]) <- 'RNA'
  test_list[[i]] <- suppressWarnings(SCTransform(test_list[[i]], verbose = TRUE))
}
for (i in 1:length(test_list)) {
  DefaultAssay(test_list[[i]]) <- 'SCT'
}
test_int_features <- SelectIntegrationFeatures(
  object.list = test_list, nfeatures = 3000)
test_list <- PrepSCTIntegration(
  object.list = test_list, anchor.features = test_int_features, 
  verbose = FALSE)
test_list <- lapply(X = test_list, FUN = function(x) {
  x <- RunPCA(x, features = test_int_features, verbose = FALSE)
})
test_anchors <- FindIntegrationAnchors(
  object.list = test_list, reduction = "rpca",
  anchor.features = test_int_features,
  scale=FALSE, normalization.method = 'SCT',
  dims = 1:30)

rm(test_list)

test_rna_integrated <- IntegrateData(anchorset = test_anchors, dims = 1:30)
test_assay <- CreateSeuratObject(
  counts = test_rna_integrated@assays$integrated,
  assay = 'RNA',
  meta.data = test_rna_integrated@meta.data,
  min.cells=0,
  min.features=0)

rm(test_rna_integrated)

# Begin Actual Modality Imputation
anchors <- FindTransferAnchors(
  reference = train_assay,
  query = test_assay,
  reference.reduction = "spca",
  dims = 1:50
)
test_assay <- TransferData(
  anchorset = anchors, 
  reference = train_assay,
  query = test_assay,
  refdata = list(
    cluster = "wsnn_res.0.2",
    celltype.l2 = 'celltype.l2',
    celltype.l3 = 'celltype.l3',
    predictedADT = "ADT")
)
test_assay <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = train_assay,
  query = test_assay, 
  new.reduction.name = "ref.spca"
)
test_assay <- ProjectUMAP(
  query = test_assay, 
  query.reduction = "ref.spca", 
  reference = train_assay, 
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

## Get corrected true counts
# Evaluate Performance
write.csv(test_assay@assays$prediction.score.celltype.l2@data, 'seurat_to_seuratlabels.csv')
write.csv(test_assay@assays$prediction.score.celltype.l3@data, 'seurat_to_seuratlabelsl3.csv')
write.csv(as.matrix(test_assay@assays$predictedADT@data), 'seurat4_protein_pbmctopbmc.csv')

protein_matrix = reference@assays$ADT@counts[, !bools]
protein_matrix = protein_matrix[, cell_indices]

protein_test = CreateSeuratObject(
  protein_matrix,
  assay = "ADT",
  meta.data = meta_data)

rm(protein_matrix)

keep_features = rownames(protein_test@assays$ADT$counts)[rowSums(protein_test@assays$ADT$counts) >= 1]
protein_test = subset(protein_test, features = keep_features)
protein_test = NormalizeData(protein_test, normalization.method = 'CLR', margin = 2) %>% ScaleData()
DimPlot(test_assay, reduction = 'ref.umap', group.by = 'predicted.celltype.l2', label.size = 3)
dim(protein_test@assays$ADT$data) #  224 85660
# subset predictions
corrs = matrix(0, nrow = nrow(protein_test@assays$ADT$data), ncol = length(unique(protein_test@meta.data$donor)))
rownames(corrs) = protein_test@assays$ADT$counts@Dimnames[1][[1]]
colnames(corrs) = unique(protein_test@meta.data$donor)
mses = corrs

keep = rep(T, nrow(corrs))

for(i in c(1:nrow(corrs))) {
  for(j in c(1:ncol(corrs))) {
    protein_name = rownames(corrs)[i]
    
    keep[i] = protein_name %in% rownames(test_assay@assays$predictedADT)
    
    if(keep[i]) {
      i_predicted = which(rownames(test_assay@assays$predictedADT) == protein_name)
      
      truth = protein_test@assays$ADT$data[i,protein_test@meta.data$donor == unique(protein_test@meta.data$donor)[j]]
      truth = truth[order(names(truth))]
      
      predicted = test_assay@assays$predictedADT$data[i_predicted,test_assay@meta.data$donor == unique(protein_test@meta.data$donor)[j]]
      predicted = predicted[order(names(predicted))]
      
      corrs[i, j] = cor(predicted, truth)
      truth = (truth - mean(truth))/sd(truth)
      predicted = (predicted - mean(predicted))/sd(predicted)
      mses[i, j] = mean((truth - predicted)^2)
    }
  }
}

corrs[is.na(corrs)] = 0
corrs = corrs[keep,]
mses = mses[keep,]
corrs = corrs[, order(colnames(corrs))]
sapply(colSums(corrs)/dim(corrs)[1], function(x) round(x, digits = 3))
mean(colSums(corrs)/dim(corrs)[1])
mses = mses[, order(colnames(mses))]
sapply(colSums(mses)/dim(mses)[1], function(x) round(x, digits = 3))
mean(colSums(mses)/dim(mses)[1])
mean(colSums(mses)/dim(mses)[1])








































































































