## This script is to generate the model with tSNE for two nonlinear clusters data
library(quollr)
library(tidyverse)

conflicted::conflicts_prefer(dplyr::filter)

source("scripts/additional_functions.R")
set.seed(20240110)

clr_choice <- "#0077A3"

## Data
training_data_two_curvy <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_data.rds") |>
  mutate(ID = row_number())

data_two_curvy <- training_data_two_curvy |>
  select(-ID) |>
  mutate(type = "data")

phate_two_curvy <- read_rds(file = "data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_phate_knn_5.rds") |>
  mutate(ID = row_number()) |>
  rename(c("emb1" = "PHATE1",
           "emb2" = "PHATE2"))

## Fit the model

num_bins_x_two_curvy <- 24

algo_obj_two_curvy <- fit_highd_model(
  highd_data = training_data_two_curvy,
  nldr_data = phate_two_curvy,
  b1 = num_bins_x_two_curvy,
  q = 0.1,
  hd_thresh = 0)

phate_two_curvy_scaled <- algo_obj_two_curvy$nldr_scaled_obj$scaled_nldr
tr_from_to_df_two_curvy <- algo_obj_two_curvy$trimesh_data
df_bin_centroids_two_curvy <- algo_obj_two_curvy$model_2d
df_bin_two_curvy <- algo_obj_two_curvy$model_highd
hex_grid <- algo_obj_two_curvy$hb_obj$hex_poly
counts_df <- algo_obj_two_curvy$hb_obj$std_cts

## Added true model

true_model_df <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_true_model.rds")
wireframe_true_model <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_true_model_connections.rds")

true_model_two_curvy <- true_model_df |>
  mutate(type = "true model")

df_bin_two_curvy <- df_bin_two_curvy |>
  select(-h) |>
  mutate(type = "model")

# Apply the scaling
df_model_data_two_curvy <- bind_rows(data_two_curvy, true_model_two_curvy, df_bin_two_curvy)

scaled_two_curvy <- scale_data_manual(df_model_data_two_curvy, "type") |>
  as_tibble()

scaled_two_curvy_data <- scaled_two_curvy |>
  filter(type == "data") |>
  select(-type)

scaled_two_curvy_data_true_model <- scaled_two_curvy |>
filter(type == "true model") |>
  select(-type)

scaled_two_curvy_data_model <- scaled_two_curvy |>
  filter(type == "model") |>
  select(-type)

## Visualize data + model

df_model_data_two_curvy_filtered <- bind_rows(df_bin_two_curvy, data_two_curvy)

langevitour::langevitour(df_model_data_two_curvy_filtered[1:(length(df_model_data_two_curvy_filtered)-1)],
                         lineFrom = tr_from_to_df_two_curvy$from,
                         lineTo = tr_from_to_df_two_curvy$to,
                         group = factor(df_model_data_two_curvy_filtered$type,
                                        c("data", "model")),
                         levelColors = c(clr_choice, "#000000"))

