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
data_two_curvy <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_data.rds") |>
  mutate(ID = row_number())

## To generate cluster column to make the split
data_two_curvy <- data_two_curvy |>
  mutate(cluster = if_else(ID <= 1000, "cluster1", "cluster2"))

## To split the data

split <- initial_split(data_two_curvy, prop = 0.7, strata = cluster)
training_data_two_curvy <- training(split)
testing_data_two_curvy <- testing(split)

write_rds(testing_data_two_curvy, file = "data/two_nonlinear/testing_two_non_linear_diff_shaped_close_clusters.rds")

################################################################################

## UMAP

n_neighbors <- 15
min_dist <- 0.1

# Create a config list with the desired parameters
umap_config <- umap.defaults
umap_config$n_neighbors <- n_neighbors      # Set the number of neighbors
umap_config$n_components <- 2    # Set the number of output dimensions (typically 2 or 3)
umap_config$min_dist <- min_dist

UMAP_fit <- umap(training_data_two_curvy[, 1:7], config = umap_config)

UMAP_data <- UMAP_fit$layout |>
  as_tibble()

names(UMAP_data) <- c("UMAP1", "UMAP2")

## Run only once
write_rds(UMAP_data, file = paste0("data/two_nonlinear/test_two_non_linear_diff_shaped_close_clusters_umap_n-neigbors_", n_neighbors, "_min-dist_", min_dist, ".rds"))

## Predict for true model

predict_UMAP_df <- predict(UMAP_fit, testing_data_two_curvy[,1:7]) |>
  as_tibble()

names(predict_UMAP_df) <- c("UMAP1", "UMAP2")

## Run only once
write_rds(predict_UMAP_df, file = "data/two_nonlinear/test_two_non_linear_diff_shaped_close_clusters_umap_predict.rds")

################################################################################

## Fit the model

num_bins_x_two_curvy <- 24

algo_obj_two_curvy <- fit_highd_model(
  highd_data = training_data_two_curvy[, 1:8],
  nldr_data = UMAP_data,
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

predict_df <- predict_emb(highd_data = testing_data_two_curvy[,1:8],
            model_highd = df_bin_two_curvy,
            model_2d = df_bin_centroids_two_curvy)

## Run only once
write_rds(predict_df, file = "data/two_nonlinear/test_two_non_linear_diff_shaped_close_clusters_highd_vis_predict.rds")

################################################################################

## First projection of test data
projection <- cbind(
  c(0.09800,0.01534,0.01887,0.00252,0.01737,-0.06895,-0.00886),
  c(-0.05248,-0.05845,0.06057,-0.00352,0.01697,-0.06938,0.01953))

highd_data <- dplyr::select(testing_data_two_curvy[,1:7], tidyselect::starts_with("x"))
projected <- as.matrix(testing_data_two_curvy[,1:7]) %*% projection
projected_df <- dplyr::mutate(dplyr::rename(tibble::as_tibble(projected,
                                                              .name_repair = "unique"),
                                            c(proj1 = "...1", proj2 = "...2")),
                              ID = dplyr::row_number())

axis_param = list(limits = 0.8,
                  axis_scaled = 5,
                  axis_pos_x = -0.95,
                  axis_pos_y = -0.95,
                  threshold = 0.042)

limits <- axis_param$limits
axis_scaled <- axis_param$axis_scaled
axis_pos_x <- axis_param$axis_pos_x
axis_pos_y <- axis_param$axis_pos_y
threshold <- axis_param$threshold
axes_obj <- gen_axes(proj = projection * axis_scaled, limits = limits,
                     axis_pos_x = axis_pos_x, axis_pos_y = axis_pos_y, axis_labels = names(highd_data),
                     threshold = threshold)
axes <- axes_obj$axes
circle <- axes_obj$circle

proj_obj1 <- list(projected_df = projected_df,
     axes = axes, circle = circle)

write_rds(proj_obj1, "data/two_nonlinear/test_two_nonlinear_proj_obj1.rds")
