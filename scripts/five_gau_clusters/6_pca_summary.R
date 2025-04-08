## This script is to check PCA results for Five Gaussian clusters data
library(tidyverse)

calculate_pca <- function(feature_dataset){
  pcaY_cal <- prcomp(feature_dataset, center = TRUE, scale = FALSE)
  PCAresults <- data.frame(pcaY_cal$x[, 1:4])
  summary_pca <- summary(pcaY_cal)
  var_explained_df <- data.frame(PC= paste0("PC",1:4),
                                 var_explained=(pcaY_cal$sdev[1:2])^2/sum((pcaY_cal$sdev[1:2])^2))
  return(list(prcomp_out = pcaY_cal,pca_components = PCAresults, summary = summary_pca, var_explained_pca  = var_explained_df))
}

data <- read_rds("data/five_gau_clusters/data_five_gau_with_clusts.rds")

pca_ref_calc <- calculate_pca(data |> dplyr::select(where(is.numeric)))
data_pca <- pca_ref_calc$pca_components |>
  dplyr::mutate(cluster = data$cluster)

## PC1 Vs PC2

ggplot(data_pca, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5) +
  theme(
    aspect.ratio = 1
  )

## PC1 Vs PC3

ggplot(data_pca, aes(x = PC1, y = PC3)) +
  geom_point(alpha = 0.5) +
  theme(
    aspect.ratio = 1
  )

## PC3 Vs PC4

ggplot(data_pca, aes(x = PC3, y = PC4)) +
  geom_point(alpha = 0.5) +
  theme(
    aspect.ratio = 1
  )

## For selected cluster

data_pca_cluster1 <- data_pca |>
  dplyr::filter(cluster == "cluster1")

ggplot(data_pca_cluster1, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5) +
  theme(
    aspect.ratio = 1
  )

## PC1 Vs PC3

ggplot(data_pca_cluster1, aes(x = PC1, y = PC3)) +
  geom_point(alpha = 0.5) +
  theme(
    aspect.ratio = 1
  )

## PC3 Vs PC4

ggplot(data_pca_cluster1, aes(x = PC3, y = PC4)) +
  geom_point(alpha = 0.5) +
  theme(
    aspect.ratio = 1
  )
