## This script is to generate the model with tSNE for two nonlinear clusters data
library(quollr)
library(tidyverse)

conflicts_prefer(dplyr::filter)

source("scripts/additional_functions.R")
set.seed(20240110)

clr_choice <- "#0077A3"

## Data
training_data_two_curvy <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_data.rds") |>
  mutate(ID = row_number())

data_two_curvy <- training_data_two_curvy |>
  select(-ID) |>
  mutate(type = "data")

tsne_two_curvy <- read_rds(file = "data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_47.rds") |>
  mutate(ID = row_number()) |>
  rename(c("emb1" = "tSNE1",
           "emb2" = "tSNE2"))

## Fit the model

num_bins_x_two_curvy <- 24

algo_obj_two_curvy <- fit_highd_model(
  highd_data = training_data_two_curvy,
  nldr_data = tsne_two_curvy,
  bin1 = num_bins_x_two_curvy,
  q = 0.1,
  benchmark_highdens = 0)

tsne_two_curvy_scaled <- algo_obj_two_curvy$nldr_obj$scaled_nldr
tr_from_to_df_two_curvy <- algo_obj_two_curvy$trimesh_data
df_bin_centroids_two_curvy <- algo_obj_two_curvy$model_2d
df_bin_two_curvy <- algo_obj_two_curvy$model_highd
hex_grid <- algo_obj_two_curvy$hb_obj$hex_poly
counts_df <- algo_obj_two_curvy$hb_obj$std_cts

write_rds(tsne_two_curvy_scaled, "data/two_nonlinear/tsne_two_curvy_scaled.rds")
write_rds(df_bin_centroids_two_curvy, "data/two_nonlinear/df_bin_centroids_two_curvy.rds")

# Code to draw illustration for notation
## hexagon binning to have regular hexagons
hb_obj_notation <- hex_binning(
  nldr_obj = algo_obj_two_curvy$nldr_obj,
  bin1 = 7,
  q = 0.1)

a1_temp <- hb_obj_notation$a1
a2_temp <- hb_obj_notation$a2
l_temp <- quad(a=3, b = 2 * a2_temp, c = -(a2_temp^2 + a1_temp^2))

## Data set with all centroids
all_centroids_df_temp <- hb_obj_notation$centroids
hex_grid_temp <- hb_obj_notation$hex_poly

hex_grid_temp40 <- hex_grid_temp |>
  filter(hex_poly_id == 40)

start_pt <- all_centroids_df_temp |>
  filter(hexID == 1)
d_rect <- tibble(x1min = 0,
                 x1max = 1,
                 x2min = 0,
                 x2max = diff(algo_obj_two_curvy$nldr_obj$lim2)/diff(algo_obj_two_curvy$nldr_obj$lim1))

# To move the rectangle to ignore the overlap with the centroids
# rect_adj <- tibble(x1 = 0.03, x2 = 0.03)
rect_adj <- tibble(x1 = -0.03, x2 = 0.03)


a1 <- tibble(x = all_centroids_df_temp$c_x[4],
             xend = all_centroids_df_temp$c_x[5],
             y = all_centroids_df_temp$c_y[21],
             yend = all_centroids_df_temp$c_y[21],
             label = expression(a[1]))
a2 <- tibble(x = all_centroids_df_temp$c_x[25],
             xend = all_centroids_df_temp$c_x[25],
             y = all_centroids_df_temp$c_y[25],
             yend = all_centroids_df_temp$c_y[33],
             label = expression(a[2]))
l <- tibble(x = hex_grid_temp40$x[2],
            xend = hex_grid_temp40$x[3],
            y = hex_grid_temp40$y[2],
            yend = hex_grid_temp40$y[3],
            label = expression(l))

write_rds(hex_grid_temp, "data/two_nonlinear/hex_grid_temp.rds")
write_rds(all_centroids_df_temp, "data/two_nonlinear/all_centroids_df_tempp.rds")
write_rds(start_pt, "data/two_nonlinear/start_pt.rds")
write_rds(d_rect, "data/two_nonlinear/d_rect.rds")
write_rds(a1, "data/two_nonlinear/a2_data.rds")
write_rds(a2, "data/two_nonlinear/a1_data.rds")
write_rds(l, "data/two_nonlinear/l_data.rds")
write_rds(rect_adj, "data/two_nonlinear/rect_adj.rds")


hex_grid_with_counts <- left_join(hex_grid, counts_df, by = c("hex_poly_id" = "hexID"))

hex_grid_nonempty <- hex_grid |>
  filter(hex_poly_id %in% df_bin_centroids_two_curvy$hexID)

bin_width <- algo_obj_two_curvy$hb_obj$a1

write_rds(hex_grid_with_counts, "data/two_nonlinear/hex_grid_with_counts.rds")
write_rds(hex_grid_nonempty, "data/two_nonlinear/hex_grid_nonempty.rds")
write_rds(tr_from_to_df_two_curvy, "data/two_nonlinear/tr_from_to_df_two_curvy.rds")


## Compute error
error_df_two_curvy_abs <- augment(
  highd_data = training_data_two_curvy,
  model_2d = df_bin_centroids_two_curvy,
  model_highd = df_bin_two_curvy)

error_df_two_curvy_abs <- error_df_two_curvy_abs |>
  mutate(sqrt_row_wise_total_error = sqrt(row_wise_total_error))

error_df_two_curvy_abs <- error_df_two_curvy_abs |>
  bind_cols(tsne_two_curvy_scaled |>
              select(-ID))

error_breaks <- error_df_two_curvy_abs$sqrt_row_wise_total_error |>
  quantile(probs = seq(0, 1, 0.05)) |>
  round(3)

error_breaks <- error_breaks[2:20]

## Add error type
breaks <- c(-Inf, error_breaks, Inf)
labels <- paste0("error", sprintf("%02d", 1:(length(breaks) - 1)))

error_df_two_curvy_abs <- error_df_two_curvy_abs |>
  mutate(error_cat = cut(sqrt_row_wise_total_error,
                         breaks = breaks,
                         labels = labels,
                         right = TRUE))

write_rds(error_df_two_curvy_abs, "data/two_nonlinear/error_df_two_curvy_abs.rds")


df_bin_two_curvy <- df_bin_two_curvy |>
  select(-hexID) |>
  mutate(type = "model")

# Apply the scaling
df_model_data_two_curvy <- bind_rows(data_two_curvy, df_bin_two_curvy)

scaled_two_curvy <- scale_data_manual(df_model_data_two_curvy, "type") |>
  as_tibble()

scaled_two_curvy_data <- scaled_two_curvy |>
  filter(type == "data") |>
  select(-type)

scaled_two_curvy_data_model <- scaled_two_curvy |>
  filter(type == "model") |>
  select(-type)

df_model_data_two_curvy_filtered <- bind_rows(df_bin_two_curvy, data_two_curvy)

langevitour::langevitour(df_model_data_two_curvy_filtered[1:(length(df_model_data_two_curvy_filtered)-1)],
                         lineFrom = tr_from_to_df_two_curvy$from,
                         lineTo = tr_from_to_df_two_curvy$to,
                         group = factor(df_model_data_two_curvy_filtered$type,
                                        c("data", "model")),
                         levelColors = c(clr_choice, "#000000"))

## Model error
langevitour::langevitour(data_two_curvy[1:(length(data_two_curvy)-1)],
                         group = factor(error_df_two_curvy_abs$error_cat,
                                        c("error01", "error02", "error03", "error04",
                                          "error05", "error06", "error07", "error08",
                                          "error09", "error10", "error11", "error12",
                                          "error13", "error14", "error15", "error16",
                                          "error17", "error18", "error19", "error20")),
                         levelColors = c("#FFFFC8", "#FFFAC0", "#FEF2B3", "#FCE9A3",
                                         "#FADE8F", "#F9D378", "#F8C65D", "#F6B938",
                                         "#F5AA00", "#F49B00", "#F28A00", "#F07800",
                                         "#ED6400", "#EA4D00", "#DE3900", "#CD2A00",
                                         "#BB1B00", "#A70C00", "#93001C", "#7D0025"))

## Model projections
## First projection
model_prj1 <- cbind(
  c(0.09800,0.01534,0.01887,0.00252,0.01737,-0.06895,-0.00886),
  c(-0.05248,-0.05845,0.06057,-0.00352,0.01697,-0.06938,0.01953))

proj_obj1 <- get_projection(projection = model_prj1,
                            proj_scale = 5,
                            highd_data = scaled_two_curvy_data,
                            model_highd = scaled_two_curvy_data_model,
                            trimesh_data = tr_from_to_df_two_curvy,
                            axis_param = list(limits = 0.8,
                                              axis_scaled = 5,
                                              axis_pos_x = -0.6,
                                              axis_pos_y = -0.6,
                                              threshold = 0.042))

# Changed the axis parameters
# axis_obj <- gen_axes(
#   proj = model_prj1 * 2,
#   limits = 0.9,
#   axis_pos_x = -0.4,
#   axis_pos_y = -0.4,
#   axis_labels = names(scaled_two_curvy_data),
#   threshold = 0.05)
#
# axes <- axis_obj$axes
# circle <- axis_obj$circle
#
# proj_obj1[["axes"]] <- axes
# proj_obj1[["circle"]] <- circle
proj_obj1[["cluster"]] <- error_df_two_curvy_abs$error_cat

write_rds(proj_obj1, "data/two_nonlinear/two_nonlinear_proj_obj1.rds")

## Second projection
model_prj2 <- cbind(
  c(0.08605,0.00286,0.01767,0.07878,0.00235,-0.03622,0.00929),
  c(0.02952,0.00022,0.07916,-0.04595,0.03757,-0.00621,-0.06809))

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

#Changed the axis parametersAdd commentMore actions
axis_obj <- gen_axes(
  proj = model_prj2 * 2,
  limits = 0.7,
  axis_pos_x = -0.35,
  axis_pos_y = -0.35,
  axis_labels = names(scaled_two_curvy_data),
  threshold = 0.02)

# axes <- axis_obj$axes
# circle <- axis_obj$circle
#
# proj_obj2[["axes"]] <- axes
# proj_obj2[["circle"]] <- circle

write_rds(proj_obj2, "data/two_nonlinear/two_nonlinear_proj_obj2.rds")


## hexbin-regular-two-curvy2

num_bins_x_two_curvy <- 15

algo_obj_two_curvy2 <- fit_highd_model(
  highd_data = training_data_two_curvy,
  nldr_data = tsne_two_curvy,
  bin1 = num_bins_x_two_curvy,
  q = 0.1,
  benchmark_highdens = 5)

tr_from_to_df_two_curvy2 <- algo_obj_two_curvy2$trimesh_data
df_bin_centroids_two_curvy2 <- algo_obj_two_curvy2$model_2d
df_bin_two_curvy2 <- algo_obj_two_curvy2$model_highd

hex_grid_two_curvy2 <- algo_obj_two_curvy2$hb_obj$hex_poly
counts_df_two_curvy2 <- algo_obj_two_curvy2$hb_obj$std_cts

hex_grid_with_counts_two_curvy2 <- left_join(hex_grid_two_curvy2, counts_df_two_curvy2, by = c("hex_poly_id" = "hexID"))

write_rds(hex_grid_with_counts_two_curvy2, "data/two_nonlinear/two_nonlinear_hex_grid_with_counts_two_curvy2.rds")
write_rds(tr_from_to_df_two_curvy2, "data/two_nonlinear/two_nonlinear_tr_from_to_df_two_curvy2.rds")


## hexbin-regular-two-curvy3

num_bins_x_two_curvy <- 48

algo_obj_two_curvy3 <- fit_highd_model(
  highd_data = training_data_two_curvy,
  nldr_data = tsne_two_curvy,
  bin1 = num_bins_x_two_curvy,
  q = 0.1,
  benchmark_highdens = 5)

tr_from_to_df_two_curvy3 <- algo_obj_two_curvy3$trimesh_data
df_bin_centroids_two_curvy3 <- algo_obj_two_curvy3$model_2d
df_bin_two_curvy3 <- algo_obj_two_curvy3$model_highd

hex_grid_two_curvy3 <- algo_obj_two_curvy3$hb_obj$hex_poly
counts_df_two_curvy3 <- algo_obj_two_curvy3$hb_obj$std_cts

hex_grid_with_counts_two_curvy3 <- left_join(hex_grid_two_curvy3, counts_df_two_curvy3, by = c("hex_poly_id" = "hexID"))

write_rds(hex_grid_with_counts_two_curvy3, "data/two_nonlinear/two_nonlinear_hex_grid_with_counts_two_curvy3.rds")
write_rds(tr_from_to_df_two_curvy3, "data/two_nonlinear/two_nonlinear_tr_from_to_df_two_curvy3.rds")

### Error
error_two_curvy_umap <- read_rds("data/two_nonlinear/error_two_non_linear_diff_shaped_close_clusters_tsne.rds")

## Find the minimum RMSE when have duplicate a1
error_two_curvy_umap <- error_two_curvy_umap |>
  group_by(a1) |>
  filter(bin1 == min(bin1)) |>
  ungroup() |>
  mutate(prop_dens = 1/(b*side_length^2))

base_line_prop_dens <- error_two_curvy_umap |>
  filter(a1 == min(a1)) |>
  pull(prop_dens)

error_two_curvy_umap <- error_two_curvy_umap |>
  mutate(prop_comp = prop_dens/base_line_prop_dens)

error_two_curvy_umap <- error_two_curvy_umap |>
  mutate(prop_bins = b_non_empty/b)

write_rds(error_two_curvy_umap, "data/two_nonlinear/error_two_curvy_umap.rds")
