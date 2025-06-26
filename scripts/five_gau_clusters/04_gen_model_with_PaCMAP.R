## This script is to generate the model with UMAP for Five Gaussian clusters data
library(quollr)
library(tidyverse)

source("scripts/additional_functions.R")
set.seed(20240110)

clr_choice <- "#0077A3"

## Data
training_data_gau <- read_rds("data/five_gau_clusters/data_five_gau.rds")

data_gau <- training_data_gau |>
  select(-ID) |>
  mutate(type = "data")


pacmap_data_gau <- read_rds("data/five_gau_clusters/pacmap_data_five_gau.rds")

## Compute hexbin parameters
num_bins_x_gau1 <- 28

algo_obj_gau1 <- fit_highd_model(
  highd_data = training_data_gau,
  nldr_data = pacmap_data_gau,
  b1 = num_bins_x_gau1,
  q = 0.1,
  benchmark_highdens = 5)

pacmap_gau_scaled <- algo_obj_gau1$nldr_obj$scaled_nldr
tr_from_to_df_gau1 <- algo_obj_gau1$trimesh_data
df_bin_centroids_gau1 <- algo_obj_gau1$model_2d
df_bin_gau1 <- algo_obj_gau1$model_highd

write_rds(pacmap_gau_scaled, "data/five_gau_clusters/pacmap_gau_scaled.rds")
write_rds(tr_from_to_df_gau1, "data/five_gau_clusters/tr_from_to_df_gau1_pacmap.rds")

df_b <- df_bin_gau1 |>
  dplyr::filter(h %in% df_bin_centroids_gau1$h) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the h order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_gau1$h, df_b$h),] |>
  dplyr::select(-h)

# Apply the scaling
df_model_data <- bind_rows(data_gau, df_b)
scaled_gau <- scale_data_manual(df_model_data, "type") |>
  as_tibble()

scaled_gau_data <- scaled_gau |>
  filter(type == "data") |>
  select(-type)

scaled_gau_data_model <- scaled_gau |>
  filter(type == "model") |>
  select(-type)

df_model_data_n <- bind_rows(df_b, data_gau)

langevitour::langevitour(df_model_data_n[1:(length(df_model_data_n)-1)],
                         lineFrom = tr_from_to_df_gau1$from,
                         lineTo = tr_from_to_df_gau1$to,
                         group = factor(df_model_data_n$type,
                                        c("data", "model")),
                         levelColors = c(clr_choice, "#000000"), pointSize = 0.7)


# Model projections

### First projection
projection <- cbind(
  c(0.10394,0.08379,-0.17868,0.32960),
  c(0.35111,-0.01137,-0.09657,-0.16018))

proj_obj1 <- get_projection(projection = projection,
                            proj_scale = 2.5,
                            highd_data = scaled_gau_data,
                            model_highd = scaled_gau_data_model,
                            trimesh_data = tr_from_to_df_gau1,
                            axis_param = list(limits = 1,
                                              axis_scaled = 2,
                                              axis_pos_x = -0.72,
                                              axis_pos_y = -0.72,
                                              threshold = 0.1))

write_rds(proj_obj1, "data/five_gau_clusters/pacmap_gau_proj_obj1.rds")

## Second projection
projection <- cbind(
  c(0.30943,-0.08452,-0.23005,-0.05284),
  c(0.22463,-0.06682,0.32106,0.02452))

proj_obj2 <- get_projection(projection = projection,
                            proj_scale = 3,
                            highd_data = scaled_gau_data,
                            model_highd = scaled_gau_data_model,
                            trimesh_data = tr_from_to_df_gau1,
                            axis_param = list(limits = 1.5,
                                              axis_scaled = 2,
                                              axis_pos_x = -1,
                                              axis_pos_y = -1,
                                              threshold = 0.05))

write_rds(proj_obj2, "data/five_gau_clusters/pacmap_gau_proj_obj2.rds")
