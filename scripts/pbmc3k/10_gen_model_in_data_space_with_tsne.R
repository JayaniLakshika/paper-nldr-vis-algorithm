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
num_bins_x_pbmc <- 13
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
  q = 0.1
)

df_bin_centroids_pbmc <- pbmc_model$df_bin_centroids
df_bin_pbmc <- pbmc_model$df_bin

## Triangulate bin centroids
tr1_object_pbmc <- tri_bin_centroids(
  df_bin_centroids_pbmc, x = "c_x", y = "c_y")
tr_from_to_df_pbmc <- gen_edges(
  tri_object = tr1_object_pbmc)

# tr_from_to_df_pbmc <- tr_from_to_df_pbmc |>
#   filter(row_number() != 360)


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

benchmark_pbmc <- 0.1

tr_from_to_df_pbmc <- tr_from_to_df_pbmc |>
  filter(!(row_number() %in% c(155)))

## Hexagonal binning to have regular hexagons
hb_obj_pbmc <- hex_binning(
  data = tsne_pbmc_scaled,
  bin1 = num_bins_x_pbmc,
  r2 = r2_pbmc,
  q = 0.1)

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

## First projection
projection <- cbind(
  c(0.038104,-0.001148,0.003905,-0.003332,-0.011747,-0.002688,0.003468,-0.014321,-0.000320),
  c(-0.000902,0.041414,0.005089,-0.000149,-0.000551,0.000239,-0.008013,-0.005799,-0.001411))

gen_proj_langevitour(
  points_df = df_exe,
  projection = projection,
  edge_df = distance_df_small_edges |> select(-distance)
)

## Second projection
projection <- cbind(
  c(0.013182,-0.022240,-0.001729,0.026946,-0.005289,0.001235,-0.007784,0.017500,-0.006949),
  c(0.027138,0.011939,0.009374,-0.004582,-0.004771,-0.027107,-0.005287,-0.002325,-0.007951))

gen_proj_langevitour(
  points_df = df_exe,
  projection = projection,
  edge_df = distance_df_small_edges |> select(-distance)
)

## Third projection
projection <- cbind(
  c(-0.003540,0.025008,0.023261,-0.012782,-0.011617,-0.006796,0.009566,-0.015022,0.001164),
  c(0.038201,-0.007751,0.005064,-0.014143,0.003255,-0.001091,0.001611,-0.003674,-0.008314))

gen_proj_langevitour(
  points_df = df_exe,
  projection = projection,
  edge_df = distance_df_small_edges |> select(-distance)
)

#### With scaled data

# Apply the scaling
scaled_pbmc_data <- scale_data_manual(training_data_pbmc |> select(-ID)) |>
  as_tibble()

df_b_pbmc <- df_bin_pbmc |>
  dplyr::filter(hb_id %in% df_bin_centroids_pbmc$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b_pbmc <- df_b_pbmc[match(df_bin_centroids_pbmc$hexID, df_b_pbmc$hb_id),] |>
  dplyr::select(-hb_id) |>
  select(-type)

# Apply the scaling
scaled_pbmc_data_model <- scale_data_manual(df_b_pbmc) |>
  as_tibble()

# Combine with the true model for visualization
df <- dplyr::bind_rows(scaled_pbmc_data_model |> mutate(type = "model"),
                       scaled_pbmc_data |> mutate(type = "data"))

## Set the maximum difference as the criteria
distance_df_small_edges_pbmc <- distance_pbmc |>
  dplyr::filter(distance < benchmark_pbmc)

distance_df_small_edges_pbmc <- distance_df_small_edges_pbmc |>
  filter(!(row_number() %in% c(155)))

# Visualize with langevitour
langevitour(df |> dplyr::select(-type),
            lineFrom = distance_df_small_edges_pbmc$from,
            lineTo = distance_df_small_edges_pbmc$to,
            group = df$type,
            pointSize = append(rep(1.5, NROW(scaled_pbmc_data_model)), rep(1, NROW(scaled_pbmc_data))),
            levelColors = c("#000000", "#33a02c"),
            lineColors = rep("#33a02c", nrow(distance_df_small_edges_pbmc)))
