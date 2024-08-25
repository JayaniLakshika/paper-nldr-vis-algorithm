## Import necessary libraries
library(quollr)
library(dplyr)
library(reader)
library(langevitour)

training_data_gau <- read_rds("data/five_gau_clusters/data_five_gau_training.rds")

umap_data_gau <- read_rds("data/five_gau_clusters/umap_data_five_gau.rds")
gau1_scaled_obj <- gen_scaled_data(
  data = umap_data_gau)
umap_gau_scaled <- gau1_scaled_obj$scaled_nldr

## Compute hexbin parameters
num_bins_x_gau1 <- 44
lim1 <- gau1_scaled_obj$lim1
lim2 <- gau1_scaled_obj$lim2
r2_gau1 <- diff(lim2)/diff(lim1)

gau1_model <- fit_highd_model(
  training_data = training_data_gau,
  emb_df = umap_gau_scaled,
  bin1 = num_bins_x_gau1,
  r2 = r2_gau1,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_gau1 <- gau1_model$df_bin_centroids
df_bin_gau1 <- gau1_model$df_bin

## Triangulate bin centroids
tr1_object_gau1 <- tri_bin_centroids(
  df_bin_centroids_gau1, x = "c_x", y = "c_y")
tr_from_to_df_gau1 <- gen_edges(
  tri_object = tr1_object_gau1)

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

## Hexagonal binning to have regular hexagons
hb_obj_gau1 <- hex_binning(
  data = umap_gau_scaled,
  bin1 = num_bins_x_gau1,
  r2 = r2_gau1)

umap_data_with_hb_id <- hb_obj_gau1$data_hb_id

df_all_gau1 <- dplyr::bind_cols(training_data_gau |> dplyr::select(-ID),
                                umap_data_with_hb_id)

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
                         group = df_exe$type, pointSize = append(rep(0, NROW(df_b)), rep(0.4, NROW(df))),
                         levelColors = c("#6a3d9a", "#33a02c"))

#### With scaled data

# data_gau_labels <- read_rds("data/five_gau_clusters/data_five_gau_with_labels.rds")
#
# training_data_gau_labels <- left_join(training_data_gau, data_gau_labels) |>
#   pull(cluster)

# Apply the scaling
scaled_gau_data <- scale_data_manual(training_data_gau |> select(-ID)) |>
  as_tibble()

df_b <- df_bin_gau1 |>
  dplyr::filter(hb_id %in% df_bin_centroids_gau1$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_gau1$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id) |>
  select(-type)

# Apply the scaling
scaled_gau_data_model <- scale_data_manual(df_b) |>
  as_tibble()

# Combine with the true model for visualization
df <- dplyr::bind_rows(scaled_gau_data_model |> mutate(type = "model"),
                       scaled_gau_data |> mutate(type = "data"))

## Set the maximum difference as the criteria
distance_df_small_edges <- distance_gau1 |>
  dplyr::filter(distance < benchmark_gau1)

# Visualize with langevitour
langevitour(df |> dplyr::select(-type),
            lineFrom = distance_df_small_edges$from,
            lineTo = distance_df_small_edges$to,
            group = df$type,
            pointSize = append(rep(1.5, NROW(scaled_gau_data_model)), rep(1, NROW(scaled_gau_data))),
            levelColors = c("#000000", "#33a02c"),
            lineColors = rep("#33a02c", nrow(distance_df_small_edges)))
