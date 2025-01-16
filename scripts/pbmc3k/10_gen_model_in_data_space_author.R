## Import necessary libraries
library(quollr)
library(dplyr)
library(readr)
library(langevitour)

## Import data
training_data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50.rds")
names(training_data_pbmc) <- paste0("pc", 1:50)

training_data_pbmc <- training_data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(training_data_pbmc))

umap_pbmc <- read_rds("data/pbmc3k/pbmc_umap_30_min_dist_0.3.rds")
pbmc_scaled_obj <- gen_scaled_data(
  data = umap_pbmc)
umap_pbmc_scaled <- pbmc_scaled_obj$scaled_nldr

## Compute hexbin parameters
num_bins_x_pbmc <- 60
lim1 <- pbmc_scaled_obj$lim1
lim2 <- pbmc_scaled_obj$lim2
r2_pbmc <- diff(lim2)/diff(lim1)

pbmc_model <- fit_highd_model(
  training_data = training_data_pbmc,
  emb_df = umap_pbmc_scaled,
  bin1 = num_bins_x_pbmc,
  r2 = r2_pbmc,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "pc"
)

df_bin_centroids_pbmc <- pbmc_model$df_bin_centroids
df_bin_pbmc <- pbmc_model$df_bin

## Triangulate bin centroids
tr1_object_pbmc <- tri_bin_centroids(
  df_bin_centroids_pbmc, x = "c_x", y = "c_y")
tr_from_to_df_pbmc <- gen_edges(
  tri_object = tr1_object_pbmc)

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
  data = umap_pbmc_scaled,
  bin1 = num_bins_x_pbmc,
  r2 = r2_pbmc)

umap_data_with_hb_id <- hb_obj_pbmc$data_hb_id

df_all_pbmc <- dplyr::bind_cols(training_data_pbmc |> dplyr::select(-ID),
                                 umap_data_with_hb_id)

### Define type column
df <- df_all_pbmc |>
  dplyr::select(tidyselect::starts_with("pc")) |>
  dplyr::mutate(type = "data") ## original dataset

df_b <- df_bin_pbmc |>
  dplyr::filter(hb_id %in% df_bin_centroids_pbmc$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_pbmc$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

df_exe <- dplyr::bind_rows(df_b, df)

## Set the maximum difference as the criteria
distance_df_small_edges_pbmc <- distance_pbmc |>
  dplyr::filter(distance < benchmark_pbmc)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges_pbmc$from,
                         lineTo = distance_df_small_edges_pbmc$to,
                         group = df_exe$type, pointSize = append(rep(0, NROW(df_b)), rep(0.5, NROW(df))),
                         levelColors = c("#6a3d9a", "#33a02c"))


#### With scaled data

data_pbmc <- training_data_pbmc |>
  select(-ID) |>
  mutate(type = "data")

df_b_pbmc <- df_bin_pbmc |>
  dplyr::filter(hb_id %in% df_bin_centroids_pbmc$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b_pbmc <- df_b_pbmc[match(df_bin_centroids_pbmc$hexID, df_b_pbmc$hb_id),] |>
  dplyr::select(-hb_id)


# Apply the scaling
df_model_data_pbmc <- bind_rows(data_pbmc, df_b_pbmc)
scaled_pbmc <- scale_data_manual(df_model_data_pbmc, "type") |>
  as_tibble()

scaled_pbmc_data <- scaled_pbmc |>
  filter(type == "data") |>
  select(-type)

scaled_pbmc_data_model <- scaled_pbmc |>
  filter(type == "model") |>
  select(-type)

# Combine with the true model for visualization
df <- dplyr::bind_rows(scaled_pbmc_data_model |> mutate(type = "model"),
                       scaled_pbmc_data |> mutate(type = "data"))


# Visualize with langevitour
langevitour(df |> dplyr::select(-type),
            lineFrom = distance_df_small_edges_pbmc$from,
            lineTo = distance_df_small_edges_pbmc$to,
            group = df$type,
            pointSize = append(rep(1.5, NROW(scaled_pbmc_data_model)), rep(1, NROW(scaled_pbmc_data))),
            levelColors = c('#a65628', "#000000"),
            lineColors = rep("#000000", nrow(distance_df_small_edges_pbmc)))



