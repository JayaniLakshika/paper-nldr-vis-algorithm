## This script is to generate the model with tSNE for MNIST digit 1 data

library(quollr)
library(tidyverse)

source("scripts/additional_functions.R")
set.seed(20240110)

## Data
training_data_mnist <- read_rds("data/mnist/mnist_10_pcs_of_digit_1.rds")
names(training_data_mnist) <- paste0("x", 1:NCOL(training_data_mnist))

training_data_mnist <- training_data_mnist |>
  mutate(ID = 1:NROW(training_data_mnist))

data_mnist <- training_data_mnist |>
  select(-ID) |>
  mutate(type = "data")

## NLDR
tsne_mnist2 <- read_rds("data/mnist/mnist_tsne89.rds")

## To fit the model
num_bins_x_mnist <- 30

algo_obj_mnist <- fit_highd_model(
  highd_data = training_data_mnist,
  nldr_data = tsne_mnist2,
  b1 = num_bins_x_mnist,
  q = 0.1,
  benchmark_highdens = 1)

tsne_minst_scaled <- algo_obj_mnist$nldr_obj$scaled_nldr
tr_from_to_df_mnist <- algo_obj_mnist$trimesh_data
df_bin_centroids_mnist <- algo_obj_mnist$model_2d
df_bin_mnist <- algo_obj_mnist$model_highd

tsne_minst_scaled <- tsne_minst_scaled |>
  mutate(cluster = if_else((emb1 <= 1) & (emb1 >= 0.78) & (emb2 <= 0.62) & (emb2 >= 0.39), "small_clust", "big_clust"))

write_rds(tsne_minst_scaled, "data/mnist/mnist_tsne_minst_scaled.rds")

## Compute error
error_df <- augment(
  highd_data = training_data_mnist,
  model_2d = df_bin_centroids_mnist,
  model_highd = df_bin_mnist)

## To join embedding
error_df <- error_df |>
  bind_cols(tsne_minst_scaled |>
              select(-ID))

## To transform the data
## ggplot(data = error_df, aes(x=row_wise_abs_error)) + geom_density() ## right skewed

error_df <- error_df |>
  mutate(sqrt_row_wise_tot_error = sqrt(row_wise_total_error)) |>
  mutate(sqrt_row_wise_tot_error = standardize(sqrt_row_wise_tot_error))

write_rds(error_df, "data/mnist/mnist_error_df.rds")

## Model
hexID_mnist <- df_bin_centroids_mnist |>
  #dplyr::filter(std_counts > 0) |>
  dplyr::pull(h)

df_b_mnist <- df_bin_mnist |>
  dplyr::filter(h %in% hexID_mnist) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the h order in df_b_with_center_data
df_b_mnist <- df_b_mnist[match(hexID_mnist, df_b_mnist$h),]

df_b_mnist_temp <- df_b_mnist

df_b_mnist <- df_b_mnist |>
  dplyr::select(-h)

# Apply the scaling
df_model_data_mnist <- bind_rows(data_mnist, df_b_mnist)
scaled_mnist <- scale_data_manual(df_model_data_mnist, "type") |>
  as_tibble()

scaled_mnist_data <- scaled_mnist |>
  filter(type == "data") |>
  select(-type)

scaled_mnist_data_model <- scaled_mnist |>
  filter(type == "model") |>
  select(-type)

## Langevitour animation

data_mnist_cluster <- data_mnist |>
  select(-type) |>
  mutate(type = tsne_minst_scaled$cluster)

df_model_data_mnist_n <- bind_rows(df_b_mnist, data_mnist_cluster)

langevitour::langevitour(df_model_data_mnist_n[1:(length(df_model_data_mnist_n)-1)],
                         lineFrom = tr_from_to_df_mnist$from,
                         lineTo = tr_from_to_df_mnist$to,
                         group = factor(df_model_data_mnist_n$type,
                                        c("big_clust", "small_clust", "model")),
                         levelColors = c("#999999",'#ff7f00', "#000000"))

## Linked plots
df_exe_mnist <- comb_all_data_model(
  highd_data = training_data_mnist,
  nldr_data = tsne_minst_scaled,
  model_highd = df_b_mnist_temp,
  model_2d = df_bin_centroids_mnist)

show_link_plots(point_data = df_exe_mnist,
                edge_data = tr_from_to_df_mnist)

## Model projections

### First projection
projection <- cbind(
  c(-0.04295,0.42361,-0.36702,0.13103,-0.19678,-0.29621,0.00471,0.09323,0.18822,-0.18484),
  c(-0.36644,0.05058,0.18200,0.21860,0.13583,-0.03511,0.35339,0.15821,-0.26945,-0.27926))

proj_obj1 <- get_projection(projection = projection,
                            proj_scale = 2.5,
                            highd_data = scaled_mnist_data,
                            model_highd = scaled_mnist_data_model,
                            trimesh_data = tr_from_to_df_mnist,
                            axis_param = list(limits = 1,
                                              axis_scaled = 1.1,
                                              axis_pos_x = -0.92,
                                              axis_pos_y = -0.92,
                                              threshold = 0.065))

proj_obj1[["cluster"]] <- tsne_minst_scaled$cluster

write_rds(proj_obj1, "data/mnist/mnist_proj_obj1.rds")

## Second projection
projection <- cbind(
  c(-0.33305,0.20223,0.06303,0.11124,0.13299,-0.30093,0.34082,0.15581,-0.08864,0.33573),
  c(0.53615,0.14969,-0.09393,-0.20053,0.07992,-0.34791,0.13057,-0.06814,-0.16737,0.03716))

proj_obj2 <- get_projection(projection = projection,
                            proj_scale = 2,
                            highd_data = scaled_mnist_data,
                            model_highd = scaled_mnist_data_model,
                            trimesh_data = tr_from_to_df_mnist,
                            axis_param = list(limits = 1,
                                              axis_scaled = 1.1,
                                              axis_pos_x = -0.83,
                                              axis_pos_y = -0.83,
                                              threshold = 0.065))

proj_obj2[["cluster"]] <- tsne_minst_scaled$cluster

write_rds(proj_obj2, "data/mnist/mnist_proj_obj2.rds")

## Images with different positions

## Data with pixel values
mnist_data <- read_rds("data/mnist/mnist_digit_1.rds")

mnist_data <- mnist_data |>
  mutate(instance = row_number()) |>
  gather(pixel, value, -Label, -instance) |>
  extract(pixel, "pixel", "(\\d+)", convert = TRUE) |>
  mutate(pixel = pixel - 2, x = pixel %% 28, y = 28 - pixel %/% 28)

img_right_top <- c(868, 1465, 1787, 2006, 4314)

pixels_gathered_within <-  mnist_data |>
  filter(instance %in% img_right_top)

write_rds(pixels_gathered_within, "data/mnist/mnist_img_right_top.rds")

img_middle <- c(130, 2091, 7795, 163, 77)

pixels_gathered_within <-  mnist_data |>
  filter(instance %in% img_middle)

write_rds(pixels_gathered_within, "data/mnist/mnist_img_middle.rds")

img_right_bottom <- c(2035, 7708, 7710, 4134, 5006)

pixels_gathered_within <-  mnist_data |>
  filter(instance %in% img_right_bottom)

write_rds(pixels_gathered_within, "data/mnist/mnist_img_right_bottom.rds")

img_error_outside <- c(999, 3650, 4916, 2964, 5622)

pixels_gathered_outside <-  mnist_data |>
  filter(instance %in% img_error_outside)

write_rds(pixels_gathered_outside, "data/mnist/mnist_img_error_outside.rds")

img_error_inside <- c(5173, 6981, 7228, 6673, 3179) #, 6673, 3179, 7363, 7650, 5241, 4636

pixels_gathered_inside <-  mnist_data |>
  filter(instance %in% img_error_inside)

write_rds(pixels_gathered_inside, "data/mnist/mnist_img_error_inside.rds")


