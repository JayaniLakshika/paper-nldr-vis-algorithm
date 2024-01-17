library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)

#################################
load("data/pbmc/pbmc_3k_festem/pbmc3k_clustering_UMAP.RData")
pbmc_umap <- umap.list$Festem@cell.embeddings |>
  tibble::as_tibble() |>
  mutate(cell_label = as.character(label.list$Festem)) |>
  rename(c("UMAP1" = "UMAP_1", "UMAP2" = "UMAP_2"))


write_rds(pbmc_umap, "data/pbmc/pbmc_3k_festem/pbmc_umap.rds")

rownames(umap.list$Festem@cell.embeddings)
plots.list[[1]]

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
load("data/pbmc/pbmc_3k_festem/pbmc3k_hvggenes.RData")
pbmc <- readRDS("data/pbmc/pbmc_3k_festem/pbmc3k_final.rds")
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

umap.list <- vector("list",length(gene.list))
label.list <- vector("list",length(gene.list))
names(label.list) <- c("Festem")
names(umap.list) <- c("Festem")
plots.list <- vector("list",length(gene.list))
for (i in 1:length(gene.list)){
  pbmc <- ScaleData(pbmc,features = gene.list[[i]])
  pbmc <- RunPCA(pbmc, verbose = FALSE,features = gene.list[[i]])
  pbmc <- FindNeighbors(object = pbmc, dims = 1:15)
  pbmc <- FindClusters(object = pbmc, resolution = 1)
  # label.list[[i]] <- pbmc@active.ident
  # umap.list[[i]] <- RunUMAP(pbmc, reduction = "pca", dims = 1:15)@reductions[["umap"]]
  # umap.tmp <- umap.for.plot(umap.list[[i]],label.list[[i]])
  #
  # class_avg <- umap.tmp %>%
  #   group_by(cluster) %>%
  #   summarise(
  #     UMAP_1 = median(UMAP_1),
  #     UMAP_2 = median(UMAP_2)
  #   )
  # plots.list[[i]] <- ggplot(umap.tmp, aes(x=UMAP_1, y=UMAP_2, color=cluster)) +
  #   geom_point(cex=0.5) + theme_bw()+theme(legend.position="none") +
  #   geom_text(aes(x=UMAP_1,y = UMAP_2,label = cluster), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 4)+
  #   labs(title = names(label.list)[i])
}

pbmc_pca <- pbmc@reductions$pca@cell.embeddings |>
  tibble::as_tibble()

write_rds(pbmc_pca, "data/pbmc/pbmc_3k_festem/pbmc_pca_50.rds")

#save(gene.list,label.list,umap.list,plots.list,file = "./results/pbmc3k_clustering_UMAP.RData")

