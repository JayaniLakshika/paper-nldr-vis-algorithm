## Import necessary libraries
library(quollr)
library(dplyr)
library(readr)
library(langevitour)

set.seed(20240110)

training_data_gau <- read_rds("data/five_gau_clusters/data_five_gau_training.rds")

tsne_data_gau <- read_rds("data/five_gau_clusters/tsne_data_five_gau_61.rds")
gau1_scaled_obj <- gen_scaled_data(
  data = tsne_data_gau)
tsne_gau_scaled <- gau1_scaled_obj$scaled_nldr

## Compute hexbin parameters
num_bins_x_gau1 <- 14
lim1 <- gau1_scaled_obj$lim1
lim2 <- gau1_scaled_obj$lim2
r2_gau1 <- diff(lim2)/diff(lim1)

gau1_model <- fit_highd_model(
  training_data = training_data_gau,
  emb_df = tsne_gau_scaled,
  bin1 = num_bins_x_gau1,
  r2 = r2_gau1,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x",
  q = 0.16
)

df_bin_centroids_gau1 <- gau1_model$df_bin_centroids
df_bin_gau1 <- gau1_model$df_bin

## Triangulate bin centroids
tr1_object_gau1 <- tri_bin_centroids(
  df_bin_centroids_gau1, x = "c_x", y = "c_y")
tr_from_to_df_gau1 <- gen_edges(
  tri_object = tr1_object_gau1)

# tr_from_to_df_gau1 <- tr_from_to_df_gau1 |>
#   filter(row_number() != 76) |>
#   filter(row_number() != 34)

tr_from_to_df_gau1 <- tr_from_to_df_gau1 |>
  filter(row_number() != 115) |>
  filter(row_number() != 139)

## Compute 2D distances
distance_gau1 <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_gau1,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_gau1 <- find_lg_benchmark(
  distance_edges = distance_gau1,
  distance_col = "distance")



trimesh_removed_gau1 <- vis_rmlg_mesh(
  distance_edges = distance_gau1,
  benchmark_value = benchmark_gau1,
  tr_coord_df = tr_from_to_df_gau1,
  distance_col = "distance")


## Hexagonal binning to have regular hexagons
hb_obj_gau1 <- hex_binning(
  data = tsne_gau_scaled,
  bin1 = num_bins_x_gau1,
  r2 = r2_gau1,
  q = 0.16)

tsne_data_with_hb_id <- hb_obj_gau1$data_hb_id

df_all_gau1 <- dplyr::bind_cols(training_data_gau |> dplyr::select(-ID),
                                  tsne_data_with_hb_id)

### Define type column
df <- df_all_gau1 |>
  dplyr::select(tidyselect::starts_with("x")) |>
  dplyr::mutate(type = "data") ## original dataset

df_b <- df_bin_gau1 |>
  dplyr::filter(hb_id %in% df_bin_centroids_gau1$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_gau1$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

df_exe <- dplyr::bind_rows(df_b, df)

## Set the maximum difference as the criteria
distance_df_small_edges <- distance_gau1 |>
  dplyr::filter(distance < benchmark_gau1)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to,
                         group = df_exe$type, pointSize = append(rep(2, NROW(df_b)), rep(1, NROW(df))),
                         levelColors = c("#6a3d9a", "#33a02c"))


## First projection
projection <- cbind(
  c(-0.00215,-0.68905,-0.04778,-0.54223),
  c(0.42558,-0.23854,-0.63659,0.35753))

gen_proj_langevitour(
  points_df = df_exe,
  projection = projection,
  edge_df = distance_df_small_edges |> select(-distance)
)
