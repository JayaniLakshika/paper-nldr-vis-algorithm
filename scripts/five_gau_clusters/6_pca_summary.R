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
  pcaRotations <- data.frame(pcaY_cal$rotation[, 1:4])
  summary_pca <- summary(pcaY_cal)
  var_explained_df <- data.frame(PC= paste0("PC",1:4),
                                 var_explained=(pcaY_cal$sdev[1:2])^2/sum((pcaY_cal$sdev[1:2])^2))
  return(list(prcomp_out = pcaY_cal,pca_components = PCAresults, summary = summary_pca, var_explained_pca  = var_explained_df, rotations = pcaRotations))
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

rotations_df <- pca_ref_calc$rotations

################## Model for PaCMAP #######################################
model_model <- read_rds("data/five_gau_clusters/pacmap_model.rds") |>
  dplyr::select(x1:x4)

projected_model <- as.matrix(model_model) %*% as.matrix(rotations_df)
projected_model <- projected_model |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::mutate(ID = dplyr::row_number())

model_wireframe <- read_rds("data/five_gau_clusters/pacmap_wireframe.rds")
model_wireframe <- model_wireframe |>
  dplyr::select(from, to)

model_wireframe <- left_join(model_wireframe, projected_model, by = c("from" = "ID"))
names(model_wireframe)[3:6] <- paste0("from_", names(model_wireframe)[3:6])

model_wireframe <- left_join(model_wireframe, projected_model, by = c("to" = "ID"))
names(model_wireframe)[7:10] <- paste0("to_", names(model_wireframe)[7:10])

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
    data = model_wireframe,
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
    data = model_wireframe,
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
    data = model_wireframe,
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

pacmap_map_df <- read_rds("data/five_gau_clusters/pacmap_model_mapping_data.rds") |>
  rename("hexID" = "ID",
         "ID" = "pts_ID") |>
  select(-hb_id) |>
  arrange(ID) |>
  mutate(cluster = data$cluster) |>
  group_by(hexID) |>
  summarize(cluster_list = list(cluster), .groups = "drop")


model_wireframe <- left_join(model_wireframe, pacmap_map_df, by = c("from" = "hexID"))
names(model_wireframe)[11] <- "from_cluster_list"

model_wireframe <- left_join(model_wireframe, pacmap_map_df, by = c("to" = "hexID"))
names(model_wireframe)[12] <- "to_cluster_list"

unlisted_model_wireframe_df <- model_wireframe |>
  unnest(cols = c(from_cluster_list)) |>
  rename(c(from_cluster = from_cluster_list)) |>
  unnest(cols = c(to_cluster_list)) |>
  rename(c(to_cluster = to_cluster_list))

data_pca_cluster1 <- data_pca |>
  dplyr::filter(cluster == "cluster1")

model_pacmap_cluster1 <- unlisted_model_wireframe_df |>
  dplyr::filter(from_cluster == "cluster1") |>
  dplyr::filter(to_cluster == "cluster1") |>
  distinct()

ggplot(data_pca_cluster1, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.1) +
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
  geom_point(alpha = 0.1) +
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
  geom_point(alpha = 0.1) +
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
