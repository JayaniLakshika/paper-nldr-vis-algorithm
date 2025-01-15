## Import necessary libraries
library(quollr)
library(dplyr)
library(readr)
library(langevitour)

## Import data
training_data_mnist <- read_rds("data/mnist/mnist_10_pcs_of_digit_1.rds")
training_data_mnist <- training_data_mnist |>
  mutate(ID = 1:NROW(training_data_mnist))

tsne_minst <- read_rds("data/mnist/mnist_tsne89.rds")

mnist_scaled_obj <- gen_scaled_data(
  data = tsne_minst)
tsne_minst_scaled <- mnist_scaled_obj$scaled_nldr

## Compute hexbin parameters
num_bins_x_mnist <- 30
lim1 <- mnist_scaled_obj$lim1
lim2 <- mnist_scaled_obj$lim2
r2_mnist <- diff(lim2)/diff(lim1)

mnist_model <- fit_highd_model(
  training_data = training_data_mnist,
  emb_df = tsne_minst_scaled,
  bin1 = num_bins_x_mnist,
  r2 = r2_mnist,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = TRUE,
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

#benchmark_mnist <- 0.1

## Hexagonal binning to have regular hexagons
hb_obj_mnist <- hex_binning(
  data = tsne_minst_scaled,
  bin1 = num_bins_x_mnist,
  r2 = r2_mnist)

tsne_data_with_hb_id <- hb_obj_mnist$data_hb_id
df_all_mnist <- dplyr::bind_cols(training_data_mnist |> dplyr::select(-ID),
                           tsne_data_with_hb_id)

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
distance_df_small_edges_mnist <- distance_mnist |>
  dplyr::filter(distance < benchmark_mnist)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges_mnist$from,
                         lineTo = distance_df_small_edges_mnist$to,
                         group = df_exe$type, pointSize = append(rep(1, NROW(df_b)), rep(0.5, NROW(df))),
                         levelColors = c("#6a3d9a", "#33a02c"))

#### With scaled data

data_mnist <- training_data_mnist |>
  select(-ID) |>
  mutate(type = "data")

# Apply the scaling
df_model_data <- bind_rows(data_mnist, df_b)
scaled_mnist <- scale_data_manual(df_model_data, "type") |>
  as_tibble()

scaled_mnist_data <- scaled_mnist |>
  filter(type == "data") |>
  select(-type)

scaled_mnist_data_model <- scaled_mnist |>
  filter(type == "model") |>
  select(-type)


df_b_mnist <- df_bin_mnist |>
  dplyr::filter(hb_id %in% df_bin_centroids_mnist$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b_mnist <- df_b_mnist[match(df_bin_centroids_mnist$hexID, df_b_mnist$hb_id),] |>
  dplyr::select(-hb_id) |>
  select(-type)

# Combine with the true model for visualization
df <- dplyr::bind_rows(scaled_mnist_data_model |> mutate(type = "model"),
                       scaled_mnist_data |> mutate(type = "data"))

## Set the maximum difference as the criteria
distance_df_small_edges_mnist <- distance_mnist |>
  dplyr::filter(distance < benchmark_mnist)

# Visualize with langevitour
langevitour(df |> dplyr::select(-type),
            lineFrom = distance_df_small_edges_mnist$from,
            lineTo = distance_df_small_edges_mnist$to,
            group = df$type,
            pointSize = append(rep(1.5, NROW(scaled_mnist_data_model)), rep(1, NROW(scaled_mnist_data))),
            levelColors = c("#000000", "#33a02c"),
            lineColors = rep("#33a02c", nrow(distance_df_small_edges_mnist)))


