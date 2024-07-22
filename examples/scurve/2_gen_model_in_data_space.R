## Import necessary libraries
library(quollr)
library(dplyr)
library(reader)
library(langevitour)

training_data_scurve <- read_rds("data/s_curve/s_curve_training.rds")

umap_scurve <- read_rds(file = "data/s_curve/s_curve_umap.rds")

scurve_scaled_obj <- gen_scaled_data(
  data = umap_scurve)

umap_scurve_scaled <- scurve_scaled_obj$scaled_nldr
lim1 <- scurve_scaled_obj$lim1
lim2 <- scurve_scaled_obj$lim2
r2 <- diff(lim2)/diff(lim1)

## Compute hexbin parameters
num_bins_x_scurve <- 12

scurve_model <- fit_highd_model(
  training_data = training_data_scurve,
  emb_df = umap_scurve_scaled,
  bin1 = num_bins_x_scurve,
  r2 = r2,
  q = 0.05,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = TRUE,
  col_start_highd = "x"
)

df_bin_centroids_scurve <- scurve_model$df_bin_centroids
df_bin_scurve <- scurve_model$df_bin

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

benchmark_scurve <- 0.2

trimesh_removed_scurve <- vis_rmlg_mesh(
  distance_edges = distance_scurve,
  benchmark_value = benchmark_scurve,
  tr_coord_df = tr_from_to_df_scurve,
  distance_col = "distance")

## Hexagonal binning to have regular hexagons
hb_obj_scurve <- hex_binning(
  data = umap_scurve_scaled,
  bin1 = num_bins_x_scurve,
  r2 = r2,
  q = 0.05)

## Data set with all possible centroids in the hexagonal grid
all_centroids_df <- hb_obj_scurve$centroids
glimpse(all_centroids_df)

## Generate all coordinates of hexagons
hex_grid <- hb_obj_scurve$hex_poly
glimpse(hex_grid)

## To obtain the standardise counts within hexbins
counts_df <- hb_obj_scurve$std_cts
df_bin_centroids <- extract_hexbin_centroids(centroids_df = all_centroids_df,
                                             counts_df = counts_df) |>
  filter(drop_empty == FALSE)

hex_grid_with_counts <- left_join(hex_grid, counts_df, by = c("hex_poly_id" = "hb_id"))

ggplot(data = hex_grid_with_counts, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = hex_poly_id, fill = std_counts)) +
  geom_text(data = all_centroids_df, aes(x = c_x, y = c_y, label = hexID)) +
  scale_fill_viridis_c(direction = -1, na.value = "#ffffff") +
  coord_fixed()

umap_data_with_hb_id <- hb_obj_scurve$data_hb_id

df_all_scurve <- dplyr::bind_cols(training_data_scurve |> dplyr::select(-ID),
                                umap_data_with_hb_id)

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

##############################

##############################

true_model <- read_rds("scurve_true_model.rds") |>
  select(-ID) |>
  mutate(type = "true model")

training_data <- read_rds("data/s_curve/s_curve_training.rds") |>
  dplyr::select(-ID) |>
  dplyr::mutate(type = "data")

# Combine with the true model for visualization
df <- dplyr::bind_rows(true_model, df_b, training_data)

connections <- read_rds("scurve_true_model_wireframe.rds")

distance_df_small_edges <- distance_df_small_edges |>
  mutate(from = from + max(connections$from) + 1,
         to = to + max(connections$to)) |>
  select(from, to)

point_connect_all <- dplyr::bind_rows(connections, distance_df_small_edges)


# Visualize with langevitour
langevitour(df |> dplyr::select(-type),
            lineFrom = point_connect_all$from,
            lineTo = point_connect_all$to,
            group = df$type,
            pointSize = append(rep(2, NROW(true_model) +  NROW(df_b)),
                               rep(1, NROW(training_data))),
            levelColors = c("#6a3d9a", "#33a02c", "#969696"),
            lineColors = append(rep("#969696", length(connections$from)),
                                rep("#000000", length(distance_df_small_edges$from))))


