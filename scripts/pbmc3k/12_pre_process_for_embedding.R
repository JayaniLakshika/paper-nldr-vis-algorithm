library(Seurat)
library(SeuratData)

obj <- LoadData("pbmcsca")

## To find specific method
#obj@meta.data$Method

## To subset cells with sepcific method and cell types
obj1 <- subset(obj, subset = CellType %in% c("Cytotoxic T cell", "CD4+ T cell", "CD14+ monocyte", "B cell") & Method == "inDrops")
obj2 <- subset(obj, subset = CellType %in% c("Cytotoxic T cell", "CD4+ T cell", "CD14+ monocyte", "B cell") & Method == "Drop-seq")
obj3 <- subset(obj, subset = CellType %in% c("Cytotoxic T cell", "CD4+ T cell", "CD14+ monocyte", "B cell") & Method == "Seq-Well")

obj1 = NormalizeData(obj1)
obj1 = ScaleData(obj1)
obj1 = FindVariableFeatures(obj1)
obj1 = RunPCA(obj1)

Seurat::ElbowPlot(obj1)

obj1 <- RunUMAP(obj1, dims = 1:50)
DimPlot(obj1, reduction = "umap")

obj1 <- RunTSNE(obj1, dims = 1:50)
DimPlot(obj1, reduction = "tsne")

K = 50
start = Sys.time()
obj1 = RunTSNE(obj1)
result = scDEED(obj1, K = K, reduction.method = 'tsne', perplexity = 30, rerun = F)
#the rerun argument saves a bit of time- it tells the function that we already have the embeddings and so we do not rerun RunTSNE. When optimizing hyperparameters, the rerun argument is T (default) because the function needs to rerun embeddings at each hyperparameter setting

end = Sys.time()
end-start
#Time difference of 4.822411 mins

dubious_cells = result$full_results$dubious_cells[result$full_results$perplexity=='40']
dubious_cells = as.numeric(strsplit(dubious_cells, ',')[[1]])
trustworthy_cells =  result$full_results$trustworthy_cells[result$full_results$perplexity=='40']
trustworthy_cells = as.numeric(strsplit(trustworthy_cells, ',')[[1]])

DimPlot(data, reduction = 'tsne', cells.highlight = list('dubious' = dubious_cells, 'trustworthy' = trustworthy_cells)) + scale_color_manual(values = c('gray', 'blue', 'red'))
