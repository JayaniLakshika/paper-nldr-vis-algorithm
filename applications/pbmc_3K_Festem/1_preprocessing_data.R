library(Seurat)

# Preprocessing -----------------------------------------------------------
# Based on Seurat's tutorial. https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# Load the PBMC dataset
# The downloaded file should be placed under the same folder.
pbmc.data <- Read10X(data.dir = "data/pbmc/pbmc_3k_festem/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc <- UpdateSeuratObject(pbmc)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


# Pre-clustering ----------------------------------------------------------
## Using top 8000 genes selected by HVGvst as prior for Festem.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 8000)
pbmc <- ScaleData(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap",label = T)
cluster.labels <- pbmc@active.ident
#save(cluster.labels,file = "./results/pbmc_preclustering.RData")

pbmc <- FindClusters(pbmc, resolution = 0.8)
cluster.labels <- pbmc@active.ident
#save(cluster.labels,file = "./results/pbmc_preclustering_10g.RData")

pbmc <- FindClusters(pbmc, resolution = 0.3)
cluster.labels <- pbmc@active.ident
#save(cluster.labels,file = "./results/pbmc_preclustering_7g.RData")

## Using top 2000 genes selected by HVGvst as labels for DEG detection methods.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap",label = T)
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
cluster.labels <- pbmc@active.ident
#save(cluster.labels,file = "./results/pbmc3k_label.RData")
#DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(pbmc, file = "data/pbmc/pbmc_3k_festem/pbmc3k.rds")
