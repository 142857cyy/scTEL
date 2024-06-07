library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)


reference = LoadH5Seurat("../Data/pbmc_multimodal.h5seurat")
gene_matrix <- reference@assays$SCT@counts 
protein_matrix <- reference@assays$ADT@counts

cell_indices = colSums(gene_matrix >= 1) >= 200
gene_matrix = gene_matrix[, cell_indices] # 20729,161748
protein_matrix = protein_matrix[, cell_indices] 
meta_data = reference@meta.data[cell_indices,]
rm(reference)


train_assay = CreateSeuratObject(gene_matrix,
                                 assay = "RNA",
                                 meta.data = meta_data)
rm(gene_matrix)
protein = CreateSeuratObject(protein_matrix,
                             assay = "ADT",
                             meta.data = meta_data)
rm(protein_matrix)
keep_features = rownames(train_assay@assays$RNA$counts)[rowSums(train_assay@assays$RNA$counts >= 1) >= 30]
train_assay = subset(train_assay,features=keep_features)
keep_features = rownames(protein@assays$ADT$counts)[rowSums(protein@assays$ADT$counts) >= 1]
protein = subset(protein,features = keep_features)
train_assay@assays$ADT = protein@assays$ADT
train_list <- SplitObject(train_assay, split.by = "donor")


for (i in 1:length(train_list)) {
  train_list[[i]] <- SCTransform(train_list[[i]], verbose = TRUE)
  DefaultAssay(train_list[[i]]) <- 'ADT'
  VariableFeatures(train_list[[i]]) <- rownames(train_list[[i]][["ADT"]])
  train_list[[i]] <- NormalizeData(train_list[[i]], normalization.method = 'CLR', margin = 2) %>% 
    ScaleData()
  
}

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
SaveH5Seurat(train_rna_integrated, "train_integrated.h5seurat", overwrite = TRUE)
rm(train_rna_integrated)

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
train_rna_integrated = LoadH5Seurat("train_integrated.h5seurat")
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

# Impute Query Dataset
gene_matrix = readMM("../Data/H1N1/gene_data.mtx")  # 32738,53201
protein_matrix = readMM("../Data/H1N1/protein_data.mtx") # 87,53201
# 53201
meta_data = read.table("../Data/H1N1/meta_data.txt", sep = ',', row.names = 1, header = T)
colnames(gene_matrix) = colnames(protein_matrix) = rownames(meta_data)
rownames(gene_matrix) = read.table("../Data/H1N1/gene_names.txt", sep = ',')$V2[-1]
rownames(protein_matrix) = read.table("../Data/H1N1/protein_names.txt", sep = ',')$V2[-1]

cell_indices = colSums(gene_matrix >= 1) >= 200
print(dim(gene_matrix)) # 32738,53201
gene_matrix = gene_matrix[, cell_indices]
protein_matrix = protein_matrix[, cell_indices]
print(dim(gene_matrix)) # 32738,53200
meta_data = meta_data[cell_indices,]

test_assay = CreateSeuratObject(
  gene_matrix,
  assay = "RNA",
  meta.data = meta_data)

rm(gene_matrix)

protein_test = CreateSeuratObject(
  protein_matrix,
  assay = "ADT",
  meta.data = meta_data)

rm(protein_matrix)

keep_features = rownames(test_assay@assays$RNA$counts)[rowSums(test_assay@assays$RNA$counts >= 1) >= 30]
test_assay = subset(test_assay, features = keep_features)

keep_features = rownames(protein_test@assays$ADT$counts)[rowSums(protein_test@assays$ADT$counts) >= 1]
protein_test = subset(protein_test, features = keep_features)

test_assay@assays$ADT = protein_test@assays$ADT #comment this

test_list <- SplitObject(test_assay, split.by = "sample")
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
# Get corrected true counts
# Evaluste Performance
protein_matrix = readMM("../Data/H1N1/protein_data.mtx")
meta_data = read.table("../Data/H1N1/meta_data.txt", sep = ',', row.names = 1, header = T)
colnames(protein_matrix) = rownames(meta_data)
rownames(protein_matrix) = read.table("../Data/H1N1/protein_names.txt", sep = ',')$V2[-1]
protein_matrix = protein_matrix[, cell_indices]
meta_data = meta_data[cell_indices, ]

protein_test = CreateSeuratObject(
  protein_matrix,
  assay = "ADT",
  meta.data = meta_data)
rm(protein_matrix)
keep_features = rownames(protein_test@assays$ADT$counts)[rowSums(protein_test@assays$ADT$counts) >= 1]
protein_test = subset(protein_test, features = keep_features)
protein_test = NormalizeData(protein_test, normalization.method = 'CLR', margin = 2) %>% ScaleData()
predicted_data = test_assay@assays$predictedADT@data
dim(protein_test@assays$ADT$data) # 87,53200

# subset predictions
corrs = matrix(0, nrow = nrow(protein_test@assays$ADT$data), ncol = length(unique(protein_test@meta.data$sample))) # 87.20
rownames(corrs) = protein_test@assays$ADT$counts@Dimnames[1][[1]]
colnames(corrs) = unique(protein_test@meta.data$sample)
mses = corrs

keep = rep(T, nrow(corrs))

for(i in c(1:nrow(corrs))) {
  for(j in c(1:ncol(corrs))) {
    protein_name = rownames(corrs)[i]
    protein_name = substr(protein_name, 1, nchar(protein_name) - 5)
    
    keep[i] = protein_name %in% rownames(predicted_data)
    
    if(keep[i]) {
      i_predicted = which(rownames(predicted_data) == protein_name)
      
      truth = protein_test@assays$ADT$data[i,protein_test@meta.data$sample == unique(protein_test@meta.data$sample)[j]]
      truth = truth[order(names(truth))]
      
      predicted = predicted_data[i_predicted,test_assay@meta.data$sample == unique(protein_test@meta.data$sample)[j]]
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
corrs = corrs[, order(colnames(corrs))] # 59,20
sapply(colSums(corrs)/dim(corrs)[1], function(x) round(x, digits = 3))
mean(colSums(corrs)/dim(corrs)[1])

test_assay@assays$prediction.score.celltype.l2 # 31features,53200cells
DimPlot(test_assay, reduction = 'ref.umap', group.by = 'predicted.celltype.l2', label.size = 3)

write.csv(test_assay@assays$prediction.score.celltype.l2@data, 'seurat_to_h1n1labels.csv')

mses = mses[, order(colnames(mses))]
sapply(colSums(mses)/dim(mses)[1], function(x) round(x, digits = 3))
mean(colSums(mses)/dim(mses)[1])

rownames(corrs) = sapply(rownames(corrs), function(x) substr(x, 1, regexpr('-PROT', x) - 1))
write.csv(corrs, "corrs_results/seurat4_pbmctoh1n1.csv")

rownames(mses) = sapply(rownames(mses), function(x) substr(x, 1, regexpr('-PROT', x) - 1))
write.csv(mses, "mse_results/seurat4_pbmctoh1n1.csv")

train_assay <- RunUMAP(train_assay, dims = 1:10)

write.csv(train_assay@reductions$umap@cell.embeddings, "pbmctoh1n1_trainumap.csv")
write.csv(test_assay@reductions$ref.umap@cell.embeddings, "pbmctoh1n1_testumap.csv")
write.csv(as.matrix(test_assay@assays$predictedADT@data[rownames(corrs),]), 'seurat_pbmctoh1n1features.csv')
































