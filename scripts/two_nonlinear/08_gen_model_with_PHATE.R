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
                         lineFrom = tr_from_to_df_two_curvy$from_reindexed,
                         lineTo = tr_from_to_df_two_curvy$to_reindexed,
                         group = factor(df_model_data_two_curvy_filtered$type,
                                        c("data", "model")),
                         levelColors = c(clr_choice, "#000000"))


## First projection
model_prj1 <- cbind(
  c(0.04463,-0.02833,0.08009,0.03515,-0.05412,0.01459,-0.04136),
  c(0.04992,-0.01104,-0.00263,0.07529,0.02804,-0.02512,0.07477))

proj_obj1 <- get_projection(projection = model_prj1,
                            proj_scale = 5,
                            highd_data = scaled_two_curvy_data,
                            model_highd = scaled_two_curvy_data_model,
                            trimesh_data = tr_from_to_df_two_curvy,
                            axis_param = list(limits = 0.5,
                                              axis_scaled = 6,
                                              axis_pos_x = -0.4,
                                              axis_pos_y = -0.4,
                                              threshold = 0.042))

two_curvy_proj_tsne_all_model1 <- plot_proj(
  proj_obj = proj_obj1,
  point_param = c(0.1, 0.5, clr_choice), # size, alpha, color
  line_param = c(0.5, 0.5, "black"), # linewidth, alpha
  plot_limits = c(-0.5, 0.5),
  axis_text_size = 2,
  is_category = FALSE) +
  interior_annotation(label = "c", cex = 1.2) +
  theme(aspect.ratio = 1)

#Changed the axis parametersAdd commentMore actions
axis_obj <- gen_axes(
  proj = model_prj1 * 5,
  limits = 0.4,
  axis_pos_x = -0.33,
  axis_pos_y = -0.33,
  axis_labels = names(scaled_two_curvy_data),
  threshold = 0.026)

axes <- axis_obj$axes
circle <- axis_obj$circle

proj_obj1[["axes"]] <- axes
proj_obj1[["circle"]] <- circle

write_rds(proj_obj1, "data/two_nonlinear/two_nonlinear_proj_obj1_phate.rds")

## Second projection
model_prj2 <- cbind(
  c(0.03817,0.01561,0.04894,-0.04416,0.05023,0.00536,-0.08186),
  c(-0.04386,-0.05266,0.01757,-0.02885,-0.01841,-0.09298,-0.02181))

proj_obj2 <- get_projection(projection = model_prj2,
                            proj_scale = 5,
                            highd_data = scaled_two_curvy_data,
                            model_highd = scaled_two_curvy_data_model,
                            trimesh_data = tr_from_to_df_two_curvy,
                            axis_param = list(limits = 0.5,
                                              axis_scaled = 6,
                                              axis_pos_x = -0.4,
                                              axis_pos_y = -0.4,
                                              threshold = 0.042))

two_curvy_proj_tsne_all_model1 <- plot_proj(
  proj_obj = proj_obj2,
  point_param = c(0.1, 0.5, clr_choice), # size, alpha, color
  line_param = c(0.5, 0.5, "black"), # linewidth, alpha
  plot_limits = c(-0.5, 0.5),
  axis_text_size = 2,
  is_category = FALSE) +
  interior_annotation(label = "c", cex = 1.2) +
  theme(aspect.ratio = 1)

#Changed the axis parametersAdd commentMore actions
# axis_obj <- gen_axes(
#   proj = model_prj2 * 2,
#   limits = 0.7,
#   axis_pos_x = -0.35,
#   axis_pos_y = -0.35,
#   axis_labels = names(scaled_two_curvy_data),
#   threshold = 0.02)

# axes <- axis_obj$axes
# circle <- axis_obj$circle
#
# proj_obj2[["axes"]] <- axes
# proj_obj2[["circle"]] <- circle

write_rds(proj_obj2, "data/two_nonlinear/two_nonlinear_proj_obj2_phate.rds")

## Third projection
model_prj3 <- cbind(
  c(0.12223,-0.00018,-0.00163,-0.00991,0.01255,0.00812,-0.00582),
  c(0.00086,0.11932,-0.02268,0.00159,-0.01516,0.01600,0.00748))

proj_obj3 <- get_projection(projection = model_prj3,
                            proj_scale = 5,
                            highd_data = scaled_two_curvy_data,
                            model_highd = scaled_two_curvy_data_model,
                            trimesh_data = tr_from_to_df_two_curvy,
                            axis_param = list(limits = 0.5,
                                              axis_scaled = 6,
                                              axis_pos_x = -0.4,
                                              axis_pos_y = -0.4,
                                              threshold = 0.042))

two_curvy_proj_tsne_all_model1 <- plot_proj(
  proj_obj = proj_obj3,
  point_param = c(0.1, 0.5, clr_choice), # size, alpha, color
  line_param = c(0.5, 0.5, "black"), # linewidth, alpha
  plot_limits = c(-0.5, 0.5),
  axis_text_size = 2,
  is_category = FALSE) +
  interior_annotation(label = "c", cex = 1.2) +
  theme(aspect.ratio = 1)

#Changed the axis parametersAdd commentMore actions
axis_obj <- gen_axes(
  proj = model_prj3 * 5,
  limits = 0.6,
  axis_pos_x = -0.5,
  axis_pos_y = -0.5,
  axis_labels = names(scaled_two_curvy_data),
  threshold = 0.01)

axes <- axis_obj$axes
circle <- axis_obj$circle

proj_obj3[["axes"]] <- axes
proj_obj3[["circle"]] <- circle

write_rds(proj_obj3, "data/two_nonlinear/two_nonlinear_proj_obj3_phate.rds")


