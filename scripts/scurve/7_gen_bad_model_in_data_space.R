## Import necessary libraries
library(quollr)
library(dplyr)
library(readr)
library(langevitour)

training_data_scurve <- read_rds("data/s_curve/s_curve_training.rds")

phate_scurve <- read_rds(file = "data/s_curve/s_curve_phate.rds")

scurve_scaled_obj <- gen_scaled_data(
  data = phate_scurve)

phate_scurve_scaled <- scurve_scaled_obj$scaled_nldr
lim1 <- scurve_scaled_obj$lim1
lim2 <- scurve_scaled_obj$lim2
r2 <- diff(lim2)/diff(lim1)

## Compute hexbin parameters
num_bins_x_scurve <- 14

scurve_model <- fit_highd_model(
  training_data = training_data_scurve,
  emb_df = phate_scurve_scaled,
  bin1 = num_bins_x_scurve,
  r2 = r2,
  q = 0.06,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = TRUE,
  col_start_highd = "x"
)

df_bin_centroids_scurve <- scurve_model$df_bin_centroids
df_bin_scurve <- scurve_model$df_bin

# benchmark_rm_lwd <- 0.3
#
# df_bin_centroids_scurve <- df_bin_centroids_scurve |>
#   filter(std_counts > benchmark_rm_lwd)
#
# df_bin_scurve <- df_bin_scurve |>
#   filter(hb_id %in% df_bin_centroids_scurve$hexID)

## Triangulate bin centroids
tr1_object_scurve <- tri_bin_centroids(
  df_bin_centroids_scurve, x = "c_x", y = "c_y")
tr_from_to_df_scurve <- gen_edges(
  tri_object = tr1_object_scurve)

## Compute 2D distances
distance_scurve <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_scurve,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_scurve <- find_lg_benchmark(
  distance_edges = distance_scurve,
  distance_col = "distance")

#benchmark_scurve <- 0.19

## Hexagonal binning to have regular hexagons
hb_obj_scurve <- hex_binning(
  data = phate_scurve_scaled,
  bin1 = num_bins_x_scurve,
  r2 = r2,
  q = 0.06)

phate_data_with_hb_id <- hb_obj_scurve$data_hb_id

df_all_scurve <- dplyr::bind_cols(training_data_scurve |> dplyr::select(-ID),
                                  phate_data_with_hb_id)

### Define type column
df <- df_all_scurve |>
  dplyr::select(tidyselect::starts_with("x")) |>
  dplyr::mutate(type = "data") ## original dataset

df_b <- df_bin_scurve |>
  dplyr::filter(hb_id %in% df_bin_centroids_scurve$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_scurve$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

df_exe <- dplyr::bind_rows(df_b, df)

## Set the maximum difference as the criteria
distance_df_small_edges <- distance_scurve |>
  dplyr::filter(distance < benchmark_scurve)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to,
                         group = df_exe$type, pointSize = append(rep(0, NROW(df_b)), rep(0.5, NROW(df))),
                         levelColors = c("#6a3d9a", "#33a02c"))

bin2 <- calc_bins_y(bin1 = num_bins_x_scurve, r2 = r2, q = 0.06)$bin2

glance(
  df_bin_centroids = df_bin_centroids_scurve,
  df_bin = df_bin_scurve,
  training_data = training_data_scurve,
  newdata = NULL,
  type_NLDR = "PHATE",
  col_start = "x") |>
  mutate(bin1 = num_bins_x_scurve,
         bin2 = bin2,
         b = bin1 * bin2,
         b_non_empty = NROW(df_bin_centroids_scurve))

