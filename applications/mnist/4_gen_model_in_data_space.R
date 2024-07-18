## Import necessary libraries
library(quollr)
library(dplyr)
library(reader)
library(langevitour)

## Import data
training_data_mnist <- read_rds("data/mnist/mnist_10_pcs_of_digit_1.rds")
training_data_mnist <- training_data_mnist |>
  mutate(ID = 1:NROW(training_data_mnist))

pacmap_minst <- read_rds("data/mnist/mnist_pacmap.rds") |>
  select(PaCMAP1, PaCMAP2, ID)
mnist_scaled_obj <- gen_scaled_data(
  data = pacmap_minst)
pacmap_minst_scaled <- mnist_scaled_obj$scaled_nldr

## Compute hexbin parameters
num_bins_x_mnist <- 22
lim1 <- mnist_scaled_obj$lim1
lim2 <- mnist_scaled_obj$lim2
r2_mnist <- diff(lim2)/diff(lim1)

mnist_model <- fit_highd_model(
  training_data = training_data_mnist,
  emb_df = pacmap_minst_scaled,
  bin1 = num_bins_x_mnist,
  r2 = r2_mnist,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "PC"
)

df_bin_centroids_mnist <- mnist_model$df_bin_centroids
df_bin_mnist <- mnist_model$df_bin

## Triangulate bin centroids
tr1_object_mnist <- tri_bin_centroids(
  df_bin_centroids_mnist, x = "c_x", y = "c_y")
tr_from_to_df_mnist <- gen_edges(
  tri_object = tr1_object_mnist)

## Compute 2D distances
distance_mnist <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_mnist,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_mnist <- find_lg_benchmark(
  distance_edges = distance_mnist,
  distance_col = "distance")

## Hexagonal binning to have regular hexagons
hb_obj_mnist <- hex_binning(
  data = pacmap_minst_scaled,
  bin1 = num_bins_x_mnist,
  r2 = r2_mnist)

pacmap_data_with_hb_id <- hb_obj_mnist$data_hb_id
df_all_mnist <- dplyr::bind_cols(training_data_mnist |> dplyr::select(-ID),
                           pacmap_data_with_hb_id)

### Define type column
df <- df_all_mnist |>
  dplyr::select(tidyselect::starts_with("PC")) |>
  dplyr::mutate(type = "data") ## original dataset

df_b <- df_bin_mnist |>
  dplyr::filter(hb_id %in% df_bin_centroids_mnist$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_mnist$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

df_exe <- dplyr::bind_rows(df_b, df)

## Set the maximum difference as the criteria
distance_df_small_edges <- distance_mnist |>
  dplyr::filter(distance < benchmark_mnist)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to,
                         group = df_exe$type, pointSize = append(rep(0, NROW(df_b)), rep(0.5, NROW(df))),
                         levelColors = c("#000000", "#33a02c"))
