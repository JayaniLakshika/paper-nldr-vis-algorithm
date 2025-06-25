## This script is to check PCA results for Five ellipsessian clusters data
library(tidyverse)
library(quollr)
library(patchwork)
library(cardinalR)

# library(dplyr)
# library(purrr) ## map function
# library(rsample)
# library(ggplot2)
# library(readr)
# library(geozoo)
# library(mvtnorm)
#
# library(Rtsne)
# library(umap)
# library(phateR)
# library(reticulate)
# library(patchwork)
#
# library(grid)


set.seed(20240110)

# use_python("~/miniforge3/envs/pcamp_env/bin/python")
# use_condaenv("pcamp_env")
#
# reticulate::source_python(here::here("scripts/function_scripts/Fit_PacMAP_code.py"))
# reticulate::source_python(here::here("scripts/function_scripts/Fit_TriMAP_code.py"))
#
# source(here::here("scripts/nldr_code.R"))

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

# Stretch factors along each dimension (larger = more stretch)
stretch_factors <- c(2, 0.5, 2, 0.5)

# Create diagonal covariance matrix
Sigma <- diag(stretch_factors^2)

df1 <- gen_gaussian(n = 500, p = 4, m = c(0, 0, 0, 0), s = Sigma)
df2 <- gen_gaussian(n = 500, p = 4, m = c(0.5 * 6 ^2, 0, 0, 0), s = Sigma)

data <- dplyr::bind_rows(df1, df2) |>
  dplyr::mutate(ID = row_number())
#langevitour::langevitour(data)

### Model for PaCMAP
pacmap_df <- read_rds(here::here("data/two_ellipsoids/pacmap_data_two_ellipsoids.rds")) |>
  dplyr::mutate(ID = row_number())

num_bins_x_ellipse <- 22

algo_obj_ellipse <- gen_nldr_vis_algo_obj(
  high_d_data = data,
  nldr_data = pacmap_df,
  num_x_bins = num_bins_x_ellipse)

pacmap_ellipse_scaled <- algo_obj_ellipse$nldr_scaled
distance_ellipse <- algo_obj_ellipse$distance_df
tr_from_to_df_ellipse <- algo_obj_ellipse$tr_from_to_df
benchmark_ellipse <- algo_obj_ellipse$benchmark
df_bin_centroids_ellipse <- algo_obj_ellipse$df_bin_centroids
df_bin_ellipse <- algo_obj_ellipse$df_bin

distance_df_small_edges_ellipse <- distance_ellipse |>
  filter(distance < benchmark_ellipse)#benchmark_ellipse

### Model
model_model <- df_bin_ellipse |>
  dplyr::select(x1:x4)

### Fit PCA
pca_ref_calc <- calculate_pca(data |> dplyr::select(-ID) |> dplyr::select(where(is.numeric)))

rotations_df <- pca_ref_calc$rotations

data_pca <- as.matrix(data |> dplyr::select(-ID) |> dplyr::select(where(is.numeric))) %*% as.matrix(rotations_df)

projected_model <- as.matrix(model_model) %*% as.matrix(rotations_df)
projected_model <- projected_model |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::mutate(ID = dplyr::row_number())

model_wireframe <- distance_df_small_edges_ellipse
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

p1 + p2 + p3 +
  plot_layout(ncol = 3)

# ### Generate embeddings (Run only once)
#
# ### tSNE
# tSNE_data_ellipse <- Fit_tSNE(data, opt_perplexity = calculate_effective_perplexity(data), with_seed = 20240110)
#
# tSNE_data_ellipse <- tSNE_data_ellipse |>
#   select(-ID)
#
# plot_tSNE_2D(tSNE_data_ellipse)
#
# write_rds(tSNE_data_ellipse, file = paste0("data/two_ellipsoids/tsne_data_two_ellipsoids_", calculate_effective_perplexity(data),".rds"))
# #write_rds(tSNE_data_ellipse, file = paste0("data/two_ellipsoids/tsne_data_two_ellipsoids_30.rds"))
#
#
# ### UMAP
#
# UMAP_fit <- umap(data, n_neighbors = 15, n_components =  2)
#
# UMAP_data_ellipse <- UMAP_fit$layout |>
#   as.data.frame()
# names(UMAP_data_ellipse)[1:(ncol(UMAP_data_ellipse))] <- paste0(rep("UMAP",(ncol(UMAP_data_ellipse))), 1:(ncol(UMAP_data_ellipse)))
#
# write_rds(UMAP_data_ellipse, file = "data/two_ellipsoids/umap_data_two_ellipsoids.rds")
#
# ### TriMAP
#
# tem_dir <- tempdir()
#
# Fit_TriMAP_data(data, tem_dir)
#
# path <- file.path(tem_dir, "df_2_without_class.csv")
# path2 <- file.path(tem_dir, "dataset_3_TriMAP_values.csv")
#
# #Fit_TriMAP(as.integer(2), as.integer(5), as.integer(4), as.integer(3), path, path2)
# Fit_TriMAP(as.integer(2), as.integer(12), as.integer(4), as.integer(3), path, path2)
#
# TriMAP_data <- read_csv(path2)
#
# write_rds(TriMAP_data, file = "data/two_ellipsoids/trimap_data_two_ellipsoids.rds")
#
#
# ### PaCMAP
#
# tem_dir <- tempdir()
#
# Fit_PacMAP_data(data, tem_dir)
#
# path <- file.path(tem_dir, "df_2_without_class.csv")
# path2 <- file.path(tem_dir, "dataset_3_PaCMAP_values.csv")
#
# # Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.9, as.integer(2), path, path2)
# Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.5, as.integer(2), path, path2)
#
# PacMAP_data <- read_csv(path2)
# write_rds(PacMAP_data, file = "data/two_ellipsoids/pacmap_data_two_ellipsoids.rds")
#
#
#
# ### Phate
#
# PHATE_data <- Fit_PHATE(data, knn = 5, with_seed = 20240110)
# PHATE_data <- PHATE_data |>
#   select(PHATE1, PHATE2)
#
# write_rds(PHATE_data, file = "data/two_ellipsoids/phate_data_two_ellipsoids.rds")
#
#
#


