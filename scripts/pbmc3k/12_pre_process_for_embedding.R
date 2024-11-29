library(Seurat)
library(SeuratData)

obj <- LoadData("pbmcsca")

## To find specific method
#obj@meta.data$Method

## To subset cells with sepcific method and cell types
subset(obj, subset = CellType %in% c("Cytotoxic T cell", "CD4+ T cell", "CD14+ monocyte", "B cell") & Method == "inDrops")
subset(obj, subset = CellType %in% c("Cytotoxic T cell", "CD4+ T cell", "CD14+ monocyte", "B cell") & Method == "Drop-seq")
subset(obj, subset = CellType %in% c("Cytotoxic T cell", "CD4+ T cell", "CD14+ monocyte", "B cell") & Method == "Seq-Well")
