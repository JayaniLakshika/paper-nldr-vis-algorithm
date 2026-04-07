## This script is to generate the prediction from UMAP and our approach
library(quollr)
library(tidyverse)
library(rsample)
library(umap) #predit for uwot not working

conflicted::conflicts_prefer(dplyr::filter)

source("scripts/additional_functions.R")
set.seed(20240110)

clr_choice <- "#0077A3"

## Data
UMAP_data <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_umap_n-neigbors_15_min-dist_0.1.rds")

UMAP_data <- UMAP_data |>
  mutate(ID = 1:NROW(UMAP_data)) |>
  mutate(cluster = if_else(ID <= 1000, "cluster1", "cluster2"))

################################################################################

## To split the data

split <- initial_split(UMAP_data, prop = 0.7, strata = cluster)
training_umap_two_curvy <- training(split)
test_umap_two_curvy <- testing(split)

write_rds(training_umap_two_curvy, file = "data/two_nonlinear/training_umap_two_non_linear_diff_shaped_close_clusters.rds")
write_rds(test_umap_two_curvy, file = "data/two_nonlinear/test_umap_two_non_linear_diff_shaped_close_clusters.rds")

## Select training and test data

data_two_curvy <- read_rds(here::here("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_data.rds"))
data_two_curvy <- data_two_curvy |>
  mutate(ID = 1:NROW(data_two_curvy))


training_data_two_curvy <- data_two_curvy |>
  filter(ID %in% training_umap_two_curvy$ID)

test_data_two_curvy <- data_two_curvy |>
  filter(ID %in% test_umap_two_curvy$ID)

write_rds(training_data_two_curvy, file = "data/two_nonlinear/training_two_non_linear_diff_shaped_close_clusters.rds")
write_rds(test_data_two_curvy, file = "data/two_nonlinear/test_two_non_linear_diff_shaped_close_clusters.rds")

################################################################################

## Fit the model

num_bins_x_two_curvy <- 13

algo_obj_two_curvy <- fit_highd_model(
  highd_data = training_data_two_curvy[, 1:8],
  nldr_data = training_umap_two_curvy[, 1:3],
  b1 = num_bins_x_two_curvy,
  q = 0.1,
  hd_thresh = 0)

umap_two_curvy_scaled <- algo_obj_two_curvy$nldr_scaled_obj$scaled_nldr
tr_from_to_df_two_curvy <- algo_obj_two_curvy$trimesh_data
df_bin_centroids_two_curvy <- algo_obj_two_curvy$model_2d
df_bin_two_curvy <- algo_obj_two_curvy$model_highd
hex_grid <- algo_obj_two_curvy$hb_obj$hex_poly
counts_df <- algo_obj_two_curvy$hb_obj$std_cts

## Predict

predict_df <- predict_emb(highd_data = test_data_two_curvy[,1:8],
                          model_highd = df_bin_two_curvy,
                          model_2d = df_bin_centroids_two_curvy)

## Run only once
write_rds(predict_df, file = "data/two_nonlinear/test_two_non_linear_diff_shaped_close_clusters_highd_vis_predict_bin_13.rds")


num_bins_x_two_curvy <- 30

algo_obj_two_curvy <- fit_highd_model(
  highd_data = training_data_two_curvy[, 1:8],
  nldr_data = training_umap_two_curvy[, 1:3],
  b1 = num_bins_x_two_curvy,
  q = 0.1,
  hd_thresh = 0)

umap_two_curvy_scaled <- algo_obj_two_curvy$nldr_scaled_obj$scaled_nldr
tr_from_to_df_two_curvy <- algo_obj_two_curvy$trimesh_data
df_bin_centroids_two_curvy <- algo_obj_two_curvy$model_2d
df_bin_two_curvy <- algo_obj_two_curvy$model_highd
hex_grid <- algo_obj_two_curvy$hb_obj$hex_poly
counts_df <- algo_obj_two_curvy$hb_obj$std_cts

## Predict

predict_df <- predict_emb(highd_data = test_data_two_curvy[,1:8],
            model_highd = df_bin_two_curvy,
            model_2d = df_bin_centroids_two_curvy)

## Run only once
write_rds(predict_df, file = "data/two_nonlinear/test_two_non_linear_diff_shaped_close_clusters_highd_vis_predict_bin_30.rds")

## Fit the model

num_bins_x_two_curvy <- 54

algo_obj_two_curvy <- fit_highd_model(
  highd_data = training_data_two_curvy[, 1:8],
  nldr_data = training_umap_two_curvy[, 1:3],
  b1 = num_bins_x_two_curvy,
  q = 0.1,
  hd_thresh = 0)

umap_two_curvy_scaled <- algo_obj_two_curvy$nldr_scaled_obj$scaled_nldr
tr_from_to_df_two_curvy <- algo_obj_two_curvy$trimesh_data
df_bin_centroids_two_curvy <- algo_obj_two_curvy$model_2d
df_bin_two_curvy <- algo_obj_two_curvy$model_highd
hex_grid <- algo_obj_two_curvy$hb_obj$hex_poly
counts_df <- algo_obj_two_curvy$hb_obj$std_cts

## Predict

predict_df <- predict_emb(highd_data = test_data_two_curvy[,1:8],
                          model_highd = df_bin_two_curvy,
                          model_2d = df_bin_centroids_two_curvy)

## Run only once
write_rds(predict_df, file = "data/two_nonlinear/test_two_non_linear_diff_shaped_close_clusters_highd_vis_predict_bin_54.rds")


################################################################################

umap_test_scaled_obj <- gen_scaled_data(nldr_data = test_umap_two_curvy[, 1:3])
write_rds(umap_test_scaled_obj$scaled_nldr, file = "data/two_nonlinear/test_umap_scaled_two_non_linear_diff_shaped_close_clusters.rds")
