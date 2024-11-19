#For scDEED method: https://github.com/JSB-UCLA/scDEED?tab=readme-ov-file
# For PBMC3k data:https://github.com/XiDsLab/Festem_paper/tree/main

library(scDEED)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)

#################################
load("data/pbmc3k/pbmc3k_clustering_UMAP.RData")
pbmc_umap <- umap.list$Festem@cell.embeddings |>
  tibble::as_tibble() |>
  mutate(cell_label = as.character(label.list$Festem)) |>
  rename(c("UMAP1" = "UMAP_1", "UMAP2" = "UMAP_2"))

#write_rds(pbmc_umap, "data/pbmc/pbmc_3k_festem/pbmc_umap.rds")

rownames(umap.list$Festem@cell.embeddings)

## Obtain gene names
rw_names <- rownames(umap.list$Festem@cell.embeddings) |> str_split("-")

name_vec <- c()
for (i in 1:length(rw_names)) {
  name <- rw_names[[i]][1]
  name_vec <- append(name_vec, name)

}

####################PCA
# library(RColorBrewer)
# my.color <- hue_pal()(13)
# names(my.color) <- 1:13
# umap.for.plot <- function(umap,cluster){
#   umap <- umap@cell.embeddings
#   umap <- as.data.frame(umap)
#   cbind(umap,cluster = cluster)
# }
load("data/pbmc3k/pbmc3k_hvggenes.RData")
pbmc <- readRDS("data/pbmc3k/pbmc3k_final.rds")
pbmc <- UpdateSeuratObject(pbmc)
# ref <- pbmc@active.ident
# names_fac <- names(pbmc@active.ident)
#
# ref <- as.character(ref)
# ref[c("CGGGCATGACCCAA","CTTGATTGATCTTC")] <- "Platelet"
# ref <- as.factor(ref)
# names(ref) <- names_fac

#pbmc <- pbmc[,!colnames(pbmc) %in% c("CGGGCATGACCCAA","CTTGATTGATCTTC")]
pbmc <- pbmc[,colnames(pbmc) %in% name_vec]

#pbmc <- pbmc[,ref!="Platelet"] ## Didn't run
gene.list <- list(EM[1:1000])
# ,
# hvgvst[1:1000],
# hvgdisp[1:1000],
# dub,
# devianceFS[1:1000],
# trendvar[1:1000]

# umap.list <- vector("list",length(gene.list))
# label.list <- vector("list",length(gene.list))
# names(label.list) <- c("Festem")
# names(umap.list) <- c("Festem")
# plots.list <- vector("list",length(gene.list))
for (i in 1:length(gene.list)){
  pbmc <- RunPCA(pbmc, verbose = FALSE,features = gene.list[[i]])
}

# valid_features <- gene.list[gene.list %in% rownames(pbmc)]
#
# pbmc <- subset(pbmc, features = valid_features)

result = scDEED(pbmc, K = 9, reduction.method = 'tsne', perplexity = 30)
result_umap = scDEED(pbmc, K = 9, reduction.method = 'umap')

head(result$num_dubious)
head(result_umap$num_dubious)
