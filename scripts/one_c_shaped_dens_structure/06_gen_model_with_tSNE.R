## This script is to generate the model with tSNE for C-shaped data
library(quollr)
library(tidyverse)

source("scripts/additional_functions.R")
set.seed(20240110)

## Data
one_c_shaped_data <- read_rds(here::here("data/one_c_shaped_dens_structure/one_c_shaped_dens_data.rds"))
one_c_shaped_data <- one_c_shaped_data |>
  mutate(ID = row_number())

## NLSR
tsne_one_c_shaped <- read_rds(file = "data/one_c_shaped_dens_structure/one_c_shaped_dens_structure_tsne_perplexity_52.rds") |>
  mutate(ID = row_number())

## To fit the model
num_bins_x_one_c_shaped <- 20

algo_obj_dens_clust <- fit_highd_model(
  highd_data = one_c_shaped_data,
  nldr_data = tsne_one_c_shaped,
  bin1 = num_bins_x_one_c_shaped,
  q = 0.1,
  benchmark_highdens = 1)

tsne_one_c_shaped_scaled <- algo_obj_dens_clust$nldr_obj$scaled_nldr
tr_from_to_df_dens_clust <- algo_obj_dens_clust$trimesh_data
df_bin_centroids_dens_clust <- algo_obj_dens_clust$model_2d
df_bin_dens_clust <- algo_obj_dens_clust$model_highd

## Compute error
error_df_one_curvy_abs <- augment(
  highd_data = one_c_shaped_data,
  model_2d = df_bin_centroids_dens_clust,
  model_highd = df_bin_dens_clust)

df_bin_dens_clust_temp <- df_bin_dens_clust

error_df_one_curvy_abs_temp <- error_df_one_curvy_abs

error_df_one_curvy_abs <- error_df_one_curvy_abs |>
  bind_cols(tsne_one_c_shaped_scaled |>
              select(-ID))

error_df_one_curvy_abs <- error_df_one_curvy_abs |>
  mutate(sqrt_row_wise_total_error = sqrt(row_wise_total_error))

# Compute density
density_data <- density(error_df_one_curvy_abs$sqrt_row_wise_total_error)
density_df <- data.frame(x = density_data$x, y = density_data$y)

# Add density values to the original dataset
error_df_one_curvy_abs <- error_df_one_curvy_abs %>%
  mutate(density = approx(density_df$x, density_df$y, xout = sqrt_row_wise_total_error)$y)

error_df_one_curvy_abs <- error_df_one_curvy_abs |>
  mutate(error_cat_n = if_else(emb1 <= 0.25 & emb2 <= 0.25, "selected", "deselected")) |> ## high_error points
  mutate(error_cat_n2 = if_else(emb2 >= 1.2, "selected", "deselected")) ## corner points

error_df_one_curvy_abs_selected11 <- error_df_one_curvy_abs |>
  dplyr::filter(error_cat_n == "selected")

write_rds(error_df_one_curvy_abs_selected11, "data/one_c_shaped_dens_structure/error_df_one_curvy_abs_selected11.rds")

error_df_one_curvy_abs_deselected11 <- error_df_one_curvy_abs |>
  dplyr::filter(error_cat_n == "deselected")

write_rds(error_df_one_curvy_abs_deselected11, "data/one_c_shaped_dens_structure/error_df_one_curvy_abs_deselected11.rds")

error_df_one_curvy_abs_selected12 <- error_df_one_curvy_abs |>
  dplyr::filter(error_cat_n2 == "selected")

write_rds(error_df_one_curvy_abs_selected12, "data/one_c_shaped_dens_structure/error_df_one_curvy_abs_selected12.rds")

error_df_one_curvy_abs_deselected12 <- error_df_one_curvy_abs |>
  dplyr::filter(error_cat_n2 == "deselected")

write_rds(error_df_one_curvy_abs_deselected12, "data/one_c_shaped_dens_structure/error_df_one_curvy_abs_deselected12.rds")

# Apply the scaling

data_c_shaped <- one_c_shaped_data |>
  select(-ID) |>
  mutate(type = "data")

df_bin_dens_clust <- df_bin_dens_clust |>
  select(-hexID) |>
  mutate(type = "model")

df_model_data <- bind_rows(data_c_shaped, df_bin_dens_clust)
scaled_c_shaped <- scale_data_manual(df_model_data, "type") |>
  as_tibble()

scaled_c_shaped_data <- scaled_c_shaped |>
  filter(type == "data") |>
  select(-type)

scaled_c_shaped_data_model <- scaled_c_shaped |>
  filter(type == "model") |>
  select(-type)

## Linked plots
df_exe_c_shaped_data <- comb_all_data_model_error(
  highd_data = one_c_shaped_data,
  nldr_data = tsne_one_c_shaped,
  model_highd = df_bin_dens_clust_temp,
  model_2d = df_bin_centroids_dens_clust,
  error_data = error_df_one_curvy_abs_temp)

show_error_link_plots(point_data = df_exe_c_shaped_data,
                      edge_data = tr_from_to_df_dens_clust)

## Model projections

## First projection
projection <- cbind(
  c(0.05096,0.15399,0.19736,0.05110),
  c(0.14608,-0.16929,0.11291,-0.07160))

proj_obj1 <- get_projection(projection = projection,
                            proj_scale = 1,
                            highd_data = scaled_c_shaped_data,
                            model_highd = scaled_c_shaped_data_model,
                            trimesh_data = tr_from_to_df_dens_clust,
                            axis_param = list(limits = 0.4,
                                              axis_scaled = 3,
                                              axis_pos_x = -0.28,
                                              axis_pos_y = -0.28,
                                              threshold = 0.02))
## high error

proj_obj1[["cluster"]] <- factor(error_df_one_curvy_abs$error_cat_n2,
                                 levels=c("deselected", "selected"))

projected_df <- proj_obj1$projected_df
model_df <- proj_obj1$model_df
axes <- proj_obj1$axes
circle <- proj_obj1$circle

projected_df <- projected_df |>
  dplyr::mutate(cluster = proj_obj1$cluster)

write_rds(projected_df, "data/one_c_shaped_dens_structure/projected_df1.rds")
write_rds(model_df, "data/one_c_shaped_dens_structure/model_df1.rds")
write_rds(axes, "data/one_c_shaped_dens_structure/axes1.rds")
write_rds(circle, "data/one_c_shaped_dens_structure/circle1.rds")

## low error

proj_obj1[["cluster"]] <- factor(error_df_one_curvy_abs$error_cat_n,
                                 levels=c("deselected", "selected"))

projected_df <- proj_obj1$projected_df
model_df <- proj_obj1$model_df
axes <- proj_obj1$axes
circle <- proj_obj1$circle

projected_df <- projected_df |>
  dplyr::mutate(cluster = proj_obj1$cluster)

write_rds(projected_df, "data/one_c_shaped_dens_structure/projected_df2.rds")
write_rds(model_df, "data/one_c_shaped_dens_structure/model_df2.rds")
write_rds(axes, "data/one_c_shaped_dens_structure/axes2.rds")
write_rds(circle, "data/one_c_shaped_dens_structure/circle2.rds")

