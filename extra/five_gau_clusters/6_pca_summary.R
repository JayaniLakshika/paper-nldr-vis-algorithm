## This script is to check PCA results for Five Gaussian clusters data
library(tidyverse)
library(quollr)
library(patchwork)

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

p1 <- data_pca |>
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

p2 <- data_pca |>
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

p3 <- data_pca |>
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

## For selected cluster############################

################## Model for PaCMAP #######################################
model_model <- read_rds("data/five_gau_clusters/pacmap_model.rds") |>
  dplyr::select(x1:x4)

model_wireframe_full <- read_rds("data/five_gau_clusters/pacmap_wireframe.rds") |>
  dplyr::select(from, to)

pacmap_map_df <- read_rds("data/five_gau_clusters/pacmap_model_mapping_data.rds") |>
  rename("hexID" = "ID",
         "ID" = "pts_ID") |>
  select(-hb_id) |>
  arrange(ID) |>
  mutate(cluster = data$cluster) |>
  group_by(hexID) |>
  summarize(cluster_list = list(cluster), .groups = "drop")

# --- START: Process for a selected cluster ---
selected_cluster <- "cluster1" # Change this to the cluster you want to analyze

# 1. Filter data for the selected cluster
data_cluster_selected <- data |>
  dplyr::filter(cluster == selected_cluster) |>
  dplyr::select(where(is.numeric))

# 2. Fit PCA on the selected cluster's data
pca_cluster_calc <- calculate_pca(data_cluster_selected)
rotations_df_cluster <- pca_cluster_calc$rotations
data_pca_cluster <- as.matrix(data_cluster_selected) %*% as.matrix(rotations_df_cluster)

# 3. Project the pacmap_model onto the PCA space of the selected cluster
projected_model_cluster <- as.matrix(model_model) %*% as.matrix(rotations_df_cluster)
projected_model_cluster <- projected_model_cluster |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::mutate(Model_ID = dplyr::row_number()) # Give the projected model a temporary ID

# 4 & 5. Filter model_wireframe and join projected model coordinates
hexID_cluster_selected <- pacmap_map_df |>
  unnest(cluster_list) |>
  filter(cluster_list == selected_cluster) |>
  distinct() |>
  arrange() |>
  pull(hexID)

model_wireframe_cluster <- model_wireframe_full %>%
  dplyr::filter(from %in% hexID_cluster_selected) %>%
  dplyr::filter(to %in% hexID_cluster_selected)

model_wireframe_cluster <- left_join(model_wireframe_cluster, projected_model_cluster, by = c("from" = "Model_ID"))
names(model_wireframe_cluster)[3:6] <- paste0("from_", names(model_wireframe_cluster)[3:6])

model_wireframe_cluster <- left_join(model_wireframe_cluster, projected_model_cluster, by = c("to" = "Model_ID"))
names(model_wireframe_cluster)[7:10] <- paste0("to_", names(model_wireframe_cluster)[7:10])

p4 <- ggplot(data_pca_cluster, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.1,
             color = clr_choice) +
  geom_segment(
    data = model_wireframe_cluster,
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

p5 <- ggplot(data_pca_cluster, aes(x = PC1, y = PC3)) +
  geom_point(alpha = 0.1,
             color = clr_choice) +
  geom_segment(
    data = model_wireframe_cluster,
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

p6 <- ggplot(data_pca_cluster, aes(x = PC3, y = PC4)) +
  geom_point(alpha = 0.1,
             color = clr_choice) +
  geom_segment(
    data = model_wireframe_cluster,
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

p1 + p2 + p3 +
  p4 + p5 + p6 +
  plot_layout(ncol = 3)

################## Model for umap #######################################
model_model <- read_rds("data/five_gau_clusters/umap_model.rds") |>
  dplyr::select(x1:x4)

projected_model <- as.matrix(model_model) %*% as.matrix(rotations_df)
projected_model <- projected_model |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::mutate(ID = dplyr::row_number())

model_wireframe <- read_rds("data/five_gau_clusters/umap_wireframe.rds")
model_wireframe <- model_wireframe |>
  dplyr::select(from, to)

model_wireframe <- left_join(model_wireframe, projected_model, by = c("from" = "ID"))
names(model_wireframe)[3:6] <- paste0("from_", names(model_wireframe)[3:6])

model_wireframe <- left_join(model_wireframe, projected_model, by = c("to" = "ID"))
names(model_wireframe)[7:10] <- paste0("to_", names(model_wireframe)[7:10])

## PC1 Vs PC2

p1 <- data_pca |>
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

p2 <- data_pca |>
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

p3 <- data_pca |>
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

################## Model for umap #######################################
model_model <- read_rds("data/five_gau_clusters/umap_model.rds") |>
  dplyr::select(x1:x4)

model_wireframe_full <- read_rds("data/five_gau_clusters/umap_wireframe.rds") |>
  dplyr::select(from, to)

umap_map_df <- read_rds("data/five_gau_clusters/umap_model_mapping_data.rds") |>
  rename("hexID" = "ID",
         "ID" = "pts_ID") |>
  select(-hb_id) |>
  arrange(ID) |>
  mutate(cluster = data$cluster) |>
  group_by(hexID) |>
  summarize(cluster_list = list(cluster), .groups = "drop")

# --- START: Process for a selected cluster ---
selected_cluster <- "cluster1" # Change this to the cluster you want to analyze

# 1. Filter data for the selected cluster
data_cluster_selected <- data |>
  dplyr::filter(cluster == selected_cluster) |>
  dplyr::select(where(is.numeric))

# 2. Fit PCA on the selected cluster's data
pca_cluster_calc <- calculate_pca(data_cluster_selected)
rotations_df_cluster <- pca_cluster_calc$rotations
data_pca_cluster <- as.matrix(data_cluster_selected) %*% as.matrix(rotations_df_cluster)

# 3. Project the umap_model onto the PCA space of the selected cluster
projected_model_cluster <- as.matrix(model_model) %*% as.matrix(rotations_df_cluster)
projected_model_cluster <- projected_model_cluster |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::mutate(Model_ID = dplyr::row_number()) # Give the projected model a temporary ID

# 4 & 5. Filter model_wireframe and join projected model coordinates
hexID_cluster_selected <- umap_map_df |>
  unnest(cluster_list) |>
  filter(cluster_list == selected_cluster) |>
  distinct() |>
  arrange() |>
  pull(hexID)

model_wireframe_cluster <- model_wireframe_full %>%
  dplyr::filter(from %in% hexID_cluster_selected) %>%
  dplyr::filter(to %in% hexID_cluster_selected)

model_wireframe_cluster <- left_join(model_wireframe_cluster, projected_model_cluster, by = c("from" = "Model_ID"))
names(model_wireframe_cluster)[3:6] <- paste0("from_", names(model_wireframe_cluster)[3:6])

model_wireframe_cluster <- left_join(model_wireframe_cluster, projected_model_cluster, by = c("to" = "Model_ID"))
names(model_wireframe_cluster)[7:10] <- paste0("to_", names(model_wireframe_cluster)[7:10])

p4 <- ggplot(data_pca_cluster, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.1,
             color = clr_choice) +
  geom_segment(
    data = model_wireframe_cluster,
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

p5 <- ggplot(data_pca_cluster, aes(x = PC1, y = PC3)) +
  geom_point(alpha = 0.1,
             color = clr_choice) +
  geom_segment(
    data = model_wireframe_cluster,
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

p6 <- ggplot(data_pca_cluster, aes(x = PC3, y = PC4)) +
  geom_point(alpha = 0.1,
             color = clr_choice) +
  geom_segment(
    data = model_wireframe_cluster,
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

p1 + p2 + p3 +
  p4 + p5 + p6 +
  plot_layout(ncol = 3)


################## Model for tsne #######################################
model_model <- read_rds("data/five_gau_clusters/tsne_model.rds") |>
  dplyr::select(x1:x4)

projected_model <- as.matrix(model_model) %*% as.matrix(rotations_df)
projected_model <- projected_model |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::mutate(ID = dplyr::row_number())

model_wireframe <- read_rds("data/five_gau_clusters/tsne_wireframe.rds")
model_wireframe <- model_wireframe |>
  dplyr::select(from, to)

model_wireframe <- left_join(model_wireframe, projected_model, by = c("from" = "ID"))
names(model_wireframe)[3:6] <- paste0("from_", names(model_wireframe)[3:6])

model_wireframe <- left_join(model_wireframe, projected_model, by = c("to" = "ID"))
names(model_wireframe)[7:10] <- paste0("to_", names(model_wireframe)[7:10])

## PC1 Vs PC2

p1 <- data_pca |>
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

p2 <- data_pca |>
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

p3 <- data_pca |>
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

## For selected cluster############################

################## Model for tsne #######################################
model_model <- read_rds("data/five_gau_clusters/tsne_model.rds") |>
  dplyr::select(x1:x4)

model_wireframe_full <- read_rds("data/five_gau_clusters/tsne_wireframe.rds") |>
  dplyr::select(from, to)

tsne_map_df <- read_rds("data/five_gau_clusters/tsne_model_mapping_data.rds") |>
  rename("hexID" = "ID",
         "ID" = "pts_ID") |>
  select(-hb_id) |>
  arrange(ID) |>
  mutate(cluster = data$cluster) |>
  group_by(hexID) |>
  summarize(cluster_list = list(cluster), .groups = "drop")

# --- START: Process for a selected cluster ---
selected_cluster <- "cluster1" # Change this to the cluster you want to analyze

# 1. Filter data for the selected cluster
data_cluster_selected <- data |>
  dplyr::filter(cluster == selected_cluster) |>
  dplyr::select(where(is.numeric))

# 2. Fit PCA on the selected cluster's data
pca_cluster_calc <- calculate_pca(data_cluster_selected)
rotations_df_cluster <- pca_cluster_calc$rotations
data_pca_cluster <- as.matrix(data_cluster_selected) %*% as.matrix(rotations_df_cluster)

# 3. Project the tsne_model onto the PCA space of the selected cluster
projected_model_cluster <- as.matrix(model_model) %*% as.matrix(rotations_df_cluster)
projected_model_cluster <- projected_model_cluster |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::mutate(Model_ID = dplyr::row_number()) # Give the projected model a temporary ID

# 4 & 5. Filter model_wireframe and join projected model coordinates
hexID_cluster_selected <- tsne_map_df |>
  unnest(cluster_list) |>
  filter(cluster_list == selected_cluster) |>
  distinct() |>
  arrange() |>
  pull(hexID)

model_wireframe_cluster <- model_wireframe_full %>%
  dplyr::filter(from %in% hexID_cluster_selected) %>%
  dplyr::filter(to %in% hexID_cluster_selected)

model_wireframe_cluster <- left_join(model_wireframe_cluster, projected_model_cluster, by = c("from" = "Model_ID"))
names(model_wireframe_cluster)[3:6] <- paste0("from_", names(model_wireframe_cluster)[3:6])

model_wireframe_cluster <- left_join(model_wireframe_cluster, projected_model_cluster, by = c("to" = "Model_ID"))
names(model_wireframe_cluster)[7:10] <- paste0("to_", names(model_wireframe_cluster)[7:10])

p4 <- ggplot(data_pca_cluster, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.1,
             color = clr_choice) +
  geom_segment(
    data = model_wireframe_cluster,
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

p5 <- ggplot(data_pca_cluster, aes(x = PC1, y = PC3)) +
  geom_point(alpha = 0.1,
             color = clr_choice) +
  geom_segment(
    data = model_wireframe_cluster,
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

p6 <- ggplot(data_pca_cluster, aes(x = PC3, y = PC4)) +
  geom_point(alpha = 0.1,
             color = clr_choice) +
  geom_segment(
    data = model_wireframe_cluster,
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

p1 + p2 + p3 +
  p4 + p5 + p6 +
  plot_layout(ncol = 3)

