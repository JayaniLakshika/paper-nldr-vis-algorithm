library(getopt)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)
library(RColorBrewer)

spec <- matrix(
  c("all_analysis", "a", 0,"logical", "Perform all analysis"
  ),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

#install.packages("https://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.0.0.tar.gz", repos = NULL, type = "source")
# Download Data -----------------------------------------------------------
data_list <- InstalledData()
if (!"ifnb"%in%data_list$Dataset){
  InstallData("ifnb")
}
library(ifnb.SeuratData)
data("ifnb")
ifnb <- UpdateSeuratObject(ifnb)
ifnb <- subset(ifnb,subset = orig.ident=="IMMUNE_CTRL")
saveRDS(ifnb,file = "data/kang/ifnb_ctrl.rds")
# data("ifnb")
# ifnb <- subset(ifnb,subset = orig.ident!="IMMUNE_CTRL")
# saveRDS(ifnb,file = "./results/ifnb_stim.rds")

# system("wget --no-check-certificate https://housekeeping.unicamp.br/Housekeeping_GenesHuman.RData")
#
# if (!is.null(opt$all_analysis)){
#   system("Rscript ./3.1_preclustering.R")
# }
# system("Rscript ./3.2_run_Festem.R")
# if (!is.null(opt$all_analysis)){
#   system("Rscript ./3.3_run_DEG_methods.R")
#   system("python 3.4_run_TN_test.py")
# }
# system("Rscript ./3.5_run_FS_methods.R")
# system("Rscript ./3.6_clustering_and_tSNE.R")
# system("Rscript ./3.7_calculate_CH_indices.R")
# if (!is.null(opt$all_analysis)){
#   system("Rscript ./3.8_construct_silver_standard.R")
# }
#
# if (!is.null(opt$all_analysis)){
#   system("Rscript ./3.10_plot_figures.R -a")
# } else{
#   system("Rscript ./3.10_plot_figures.R")
# }

#################################PCA data


my.color <- hue_pal()(17)
names(my.color) <- 1:17
tsne.for.plot <- function(tsne,cluster){
  tsne <- tsne@cell.embeddings
  tsne <- as.data.frame(tsne)
  cbind(tsne,cluster = cluster)
}
load("data/kang/ifnb_ctrl_hvggenes.RData")
ifnb <- readRDS("data/kang/ifnb_ctrl.rds")
ifnb <- NormalizeData(ifnb)
ifnb <- ifnb[,!ifnb@meta.data$seurat_annotations%in%c("Eryth","Mk")]
gene.list <- list(EM[1:2500])

# ,
# hvgvst[1:2500],
# hvgdisp[1:2500],
# dub,
# devianceFS[1:2500],
# trendvar[1:2500]

tsne.list <- vector("list",length(gene.list))
label.list <- vector("list",length(gene.list))
names(label.list) <- c("Festem")
names(tsne.list) <- c("Festem")
plots.list <- vector("list",length(gene.list))
for (i in 1:length(gene.list)){
  ifnb <- ScaleData(ifnb,features = gene.list[[i]])
  ifnb <- RunPCA(ifnb, verbose = FALSE,features = gene.list[[i]])
  ifnb <- FindNeighbors(object = ifnb, dims = 1:min(25,length(gene.list[[i]])-1))
  ifnb <- FindClusters(object = ifnb, resolution = 1.5)
  # label.list[[i]] <- ifnb@active.ident
  # tsne.list[[i]] <- RunTSNE(ifnb, reduction = "pca", dims = 1:min(25,length(gene.list[[i]])-1))@reductions[["tsne"]]
  # tsne.tmp <- tsne.for.plot(tsne.list[[i]],label.list[[i]])
  #
  # class_avg <- tsne.tmp %>%
  #   group_by(cluster) %>%
  #   summarise(
  #     tSNE_1 = median(tSNE_1),
  #     tSNE_2 = median(tSNE_2)
  #   )
  # plots.list[[i]] <- ggplot(tsne.tmp, aes(x=tSNE_1, y=tSNE_2, color=cluster)) +
  #   geom_point(cex=0.5) + theme_bw()+theme(legend.position="none") +
  #   geom_text(aes(x=tSNE_1,y = tSNE_2,label = cluster), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 4)+
  #   labs(title = names(label.list)[i])
}

#save(gene.list,label.list,tsne.list,plots.list,file = "./results/ifnb_ctrl_clustering_tSNE.RData")

ifnb_pca <- ifnb@reductions$pca@cell.embeddings |>
  tibble::as_tibble()

write_rds(ifnb_pca, "data/kang/kang_pca_25.rds")


#################################tSNE data
load("data/kang/ifnb_ctrl_clustering_tSNE.RData")
kang_tsne <- tsne.list$Festem@cell.embeddings |>
  tibble::as_tibble() |>
  mutate(cell_label = as.character(label.list$Festem)) |>
  rename(c("tSNE1" = "tSNE_1", "tSNE2" = "tSNE_2"))


write_rds(kang_tsne, "data/kang/kang_tsne.rds")


plots.list[[1]]
rownames(tsne.list$Festem@cell.embeddings)
