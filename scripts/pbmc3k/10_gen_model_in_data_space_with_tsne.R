## Import necessary libraries
library(quollr)
library(dplyr)
library(reader)
library(langevitour)

## Import data
training_data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50.rds")
training_data_pbmc <- training_data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(training_data_pbmc))

tsne_pbmc <- read_rds("data/pbmc3k/pbmc_tsne_30.rds")
pbmc_scaled_obj <- gen_scaled_data(
  data = tsne_pbmc)
tsne_pbmc_scaled <- pbmc_scaled_obj$scaled_nldr

## Compute hexbin parameters
num_bins_x_pbmc <- 15
lim1 <- pbmc_scaled_obj$lim1
lim2 <- pbmc_scaled_obj$lim2
r2_pbmc <- diff(lim2)/diff(lim1)

pbmc_model <- fit_highd_model(
  training_data = training_data_pbmc,
  emb_df = tsne_pbmc_scaled,
  bin1 = num_bins_x_pbmc,
  r2 = r2_pbmc,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "PC",
  q = 0.05
)

df_bin_centroids_pbmc <- pbmc_model$df_bin_centroids
df_bin_pbmc <- pbmc_model$df_bin

## Triangulate bin centroids
tr1_object_pbmc <- tri_bin_centroids(
  df_bin_centroids_pbmc, x = "c_x", y = "c_y")
tr_from_to_df_pbmc <- gen_edges(
  tri_object = tr1_object_pbmc)

tr_from_to_df_pbmc <- tr_from_to_df_pbmc |>
  filter(row_number() != 360)

## Compute 2D distances
distance_pbmc <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_pbmc,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_pbmc <- find_lg_benchmark(
  distance_edges = distance_pbmc,
  distance_col = "distance")

## Hexagonal binning to have regular hexagons
hb_obj_pbmc <- hex_binning(
  data = tsne_pbmc_scaled,
  bin1 = num_bins_x_pbmc,
  r2 = r2_pbmc,
  q = 0.05)

tsne_data_with_hb_id <- hb_obj_pbmc$data_hb_id

df_all_pbmc <- dplyr::bind_cols(training_data_pbmc |> dplyr::select(-ID),
                                tsne_data_with_hb_id)

### Define type column
df <- df_all_pbmc |>
  dplyr::select(tidyselect::starts_with("PC")) |>
  dplyr::mutate(type = "data") ## original dataset

df_b <- df_bin_pbmc |>
  dplyr::filter(hb_id %in% df_bin_centroids_pbmc$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_pbmc$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

df_exe <- dplyr::bind_rows(df_b, df)

## Set the maximum difference as the criteria
distance_df_small_edges <- distance_pbmc |>
  dplyr::filter(distance < benchmark_pbmc)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to,
                         group = df_exe$type, pointSize = append(rep(0, NROW(df_b)), rep(0.5, NROW(df))),
                         levelColors = c("#6a3d9a", "#33a02c"))
