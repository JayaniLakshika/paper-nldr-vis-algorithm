## This script is to check PCA results for Five Gaussian clusters data
library(tidyverse)
library(quollr)
library(colorspace)

source("scripts/additional_functions.R")
set.seed(20240110)

clr_choice <- "#0077A3"

calculate_pca <- function(feature_dataset){
  pcaY_cal <- prcomp(feature_dataset, center = TRUE, scale = FALSE)
  PCAresults <- data.frame(pcaY_cal$x[, 1:4])
  summary_pca <- summary(pcaY_cal)
  var_explained_df <- data.frame(PC= paste0("PC",1:4),
                                 var_explained=(pcaY_cal$sdev[1:2])^2/sum((pcaY_cal$sdev[1:2])^2))
  return(list(prcomp_out = pcaY_cal,pca_components = PCAresults, summary = summary_pca, var_explained_pca  = var_explained_df))
}

data <- read_rds("data/five_gau_clusters/data_five_gau_with_clusts.rds")
data_n <- data
data <- data |>
  dplyr::select(-ID)

data_gau <- data |>
  select(-cluster) |>
  mutate(type = "data")

pca_ref_calc <- calculate_pca(data |> dplyr::select(where(is.numeric)))
data_pca <- pca_ref_calc$pca_components |>
  dplyr::mutate(ID = row_number()) |>
  dplyr::mutate(cluster = data$cluster)

## Model for PaCMAP
model_pacmap <- read_rds("data/five_gau_clusters/pacmap_model.rds")
model_pacmap <- model_pacmap |>
  dplyr::select(from, to)

model_pacmap <- left_join(model_pacmap, data_pca, by = c("from" = "ID"))
names(model_pacmap)[3:7] <- paste0("from_", names(model_pacmap)[3:7])

model_pacmap <- left_join(model_pacmap, data_pca, by = c("to" = "ID"))
names(model_pacmap)[8:12] <- paste0("to_", names(model_pacmap)[8:12])

## PC1 Vs PC2

data_pca |>
  ggplot(
    aes(
      x = PC1,
      y = PC2)) +
  geom_point(
    #size = 0.5,
    alpha = 0.05,
    color = clr_choice) +
  geom_segment(
    data = model_pacmap,
    aes(
      x = from_PC1,
      y = from_PC2,
      xend = to_PC1,
      yend = to_PC2),
    color = "#000000",
    #alpha = 0.4,
    linewidth = 0.5)  +
  theme(
    aspect.ratio = 1
  )

## PC1 Vs PC3

data_pca |>
  ggplot(
    aes(
      x = PC1,
      y = PC3)) +
  geom_point(
    #size = 0.5,
    alpha = 0.05,
    color = clr_choice) +
  geom_segment(
    data = model_pacmap,
    aes(
      x = from_PC1,
      y = from_PC3,
      xend = to_PC1,
      yend = to_PC3),
    color = "#000000",
    #alpha = 0.4,
    linewidth = 0.5)  +
  theme(
    aspect.ratio = 1
  )

## PC3 Vs PC4

data_pca |>
  ggplot(
    aes(
      x = PC3,
      y = PC4)) +
  geom_point(
    #size = 0.5,
    alpha = 0.05,
    color = clr_choice) +
  geom_segment(
    data = model_pacmap,
    aes(
      x = from_PC3,
      y = from_PC4,
      xend = to_PC3,
      yend = to_PC4),
    color = "#000000",
    #alpha = 0.4,
    linewidth = 0.5)  +
  theme(
    aspect.ratio = 1
  )

## For selected cluster

data_pca_cluster1 <- data_pca |>
  dplyr::filter(cluster == "cluster2")

model_pacmap_cluster1 <- model_pacmap |>
  dplyr::filter(from_cluster == "cluster2") |>
  dplyr::filter(to_cluster == "cluster2")

ggplot(data_pca_cluster1, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5) +
  geom_segment(
    data = model_pacmap_cluster1,
    aes(
      x = from_PC1,
      y = from_PC2,
      xend = to_PC1,
      yend = to_PC2),
    color = "#000000",
    #alpha = 0.4,
    linewidth = 0.5) +
  theme(
    aspect.ratio = 1
  )

## PC1 Vs PC3

ggplot(data_pca_cluster1, aes(x = PC1, y = PC3)) +
  geom_point(alpha = 0.5) +
  geom_segment(
    data = model_pacmap_cluster1,
    aes(
      x = from_PC1,
      y = from_PC3,
      xend = to_PC1,
      yend = to_PC3),
    color = "#000000",
    #alpha = 0.4,
    linewidth = 0.5) +
theme(
  aspect.ratio = 1
)

## PC3 Vs PC4

ggplot(data_pca_cluster1, aes(x = PC3, y = PC4)) +
  geom_point(alpha = 0.5) +
  geom_segment(
    data = model_pacmap_cluster1,
    aes(
      x = from_PC3,
      y = from_PC4,
      xend = to_PC3,
      yend = to_PC4),
    color = "#000000",
    #alpha = 0.4,
    linewidth = 0.5) +
theme(
  aspect.ratio = 1
)
