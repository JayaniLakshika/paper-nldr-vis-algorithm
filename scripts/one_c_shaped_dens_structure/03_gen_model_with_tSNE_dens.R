## Import necessary libraries
library(quollr)
library(dplyr)
library(readr)
library(langevitour)
library(ggplot2)

set.seed(20240110)

## Import data
training_data_one_nonlinear_structure <- read_rds("data/one_c_shaped_dens_structure/one_c_shaped_dens_data.rds")

training_data_one_nonlinear_structure <- training_data_one_nonlinear_structure |>
  mutate(ID = 1:NROW(training_data_one_nonlinear_structure))

tSNE_one_nonlinear_structure <- read_rds("data/one_c_shaped_dens_structure/one_c_shaped_dens_structure_tsne_perplexity_52.rds")

one_nonlinear_structure_scaled_obj <- gen_scaled_data(
  data = tSNE_one_nonlinear_structure)
tSNE_one_nonlinear_structure_scaled <- one_nonlinear_structure_scaled_obj$scaled_nldr |>
  mutate(ID = 1:NROW(tSNE_one_nonlinear_structure))

## Compute hexbin parameters
num_bins_x_one_nonlinear_structure <- 20
lim1 <- one_nonlinear_structure_scaled_obj$lim1
lim2 <- one_nonlinear_structure_scaled_obj$lim2
r2_one_nonlinear_structure <- diff(lim2)/diff(lim1)

one_nonlinear_structure_model <- fit_highd_model(
  training_data = training_data_one_nonlinear_structure,
  emb_df = tSNE_one_nonlinear_structure_scaled,
  bin1 = num_bins_x_one_nonlinear_structure,
  r2 = r2_one_nonlinear_structure,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_one_nonlinear_structure <- one_nonlinear_structure_model$df_bin_centroids
df_bin_one_nonlinear_structure <- one_nonlinear_structure_model$df_bin

## Triangulate bin centroids
tr1_object_one_nonlinear_structure <- tri_bin_centroids(
  df_bin_centroids_one_nonlinear_structure, x = "c_x", y = "c_y")
tr_from_to_df_one_nonlinear_structure <- gen_edges(
  tri_object = tr1_object_one_nonlinear_structure)

## Compute 2D distances
distance_one_nonlinear_structure <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_one_nonlinear_structure,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_one_nonlinear_structure <- find_lg_benchmark(
  distance_edges = distance_one_nonlinear_structure,
  distance_col = "distance")

tr_df <- distinct(tibble::tibble(
  x = c(tr_from_to_df_one_nonlinear_structure[["x_from"]], tr_from_to_df_one_nonlinear_structure[["x_to"]]),
  y = c(tr_from_to_df_one_nonlinear_structure[["y_from"]], tr_from_to_df_one_nonlinear_structure[["y_to"]])))

distance_df_small_edges_one_nonlinear_structure <- distance_one_nonlinear_structure |>
  filter(distance < benchmark_one_nonlinear_structure)

tr_from_to_df_one_nonlinear_structure <- inner_join(
  tr_from_to_df_one_nonlinear_structure, distance_df_small_edges_one_nonlinear_structure,
  by = c("from", "to"))

trimesh_removed_one_nonlinear_structure <- ggplot() +
  geom_segment(data = tr_from_to_df_one_nonlinear_structure,
               aes(
                 x = x_from,
                 y = y_from,
                 xend = x_to,
                 yend = y_to),
               colour = "#33a02c",
               linewidth = 1) +
  geom_point(data = tSNE_one_nonlinear_structure_scaled,
             aes(
               x = tSNE1,
               y = tSNE2
             ),
             alpha=0.1) +
  theme(aspect.ratio = 1)

trimesh_removed_one_nonlinear_structure


## Hexagonal binning to have regular hexagons
hb_obj_one_nonlinear_structure <- hex_binning(
  data = tSNE_one_nonlinear_structure_scaled,
  bin1 = num_bins_x_one_nonlinear_structure,
  r2 = r2_one_nonlinear_structure)

tSNE_data_with_hb_id <- hb_obj_one_nonlinear_structure$data_hb_id
df_all_one_nonlinear_structure <- dplyr::bind_cols(training_data_one_nonlinear_structure |> dplyr::select(-ID),
                                                  tSNE_data_with_hb_id)

### Define type column
df <- df_all_one_nonlinear_structure |>
  dplyr::select(tidyselect::starts_with("x")) |>
  dplyr::mutate(type = "data") ## original dataset

df_b <- df_bin_one_nonlinear_structure |>
  dplyr::filter(hb_id %in% df_bin_centroids_one_nonlinear_structure$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_one_nonlinear_structure$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

df_exe <- dplyr::bind_rows(df_b, df)

## Set the maximum difference as the criteria
# distance_df_small_edges_one_nonlinear_structure <- distance_one_nonlinear_structure |>
#   dplyr::filter(distance < benchmark_one_nonlinear_structure)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges_one_nonlinear_structure$from,
                         lineTo = distance_df_small_edges_one_nonlinear_structure$to,
                         group = df_exe$type, pointSize = append(rep(1, NROW(df_b)), rep(0.5, NROW(df))),
                         levelColors = c("#6a3d9a", "#33a02c"))

#### With scaled data

data_one_nonlinear_structure <- training_data_one_nonlinear_structure |>
  select(-ID) |>
  mutate(type = "data")

# Apply the scaling
df_model_data <- bind_rows(data_one_nonlinear_structure, df_b)
scaled_one_nonlinear_structure <- scale_data_manual(df_model_data, "type") |>
  as_tibble()

scaled_one_nonlinear_structure_data <- scaled_one_nonlinear_structure |>
  filter(type == "data") |>
  select(-type)

scaled_one_nonlinear_structure_data_model <- scaled_one_nonlinear_structure |>
  filter(type == "model") |>
  select(-type)


df_b_one_nonlinear_structure <- df_bin_one_nonlinear_structure |>
  dplyr::filter(hb_id %in% df_bin_centroids_one_nonlinear_structure$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b_one_nonlinear_structure <- df_b_one_nonlinear_structure[match(df_bin_centroids_one_nonlinear_structure$hexID, df_b_one_nonlinear_structure$hb_id),] |>
  dplyr::select(-hb_id) |>
  select(-type)

# Combine with the true model for visualization
df <- dplyr::bind_rows(scaled_one_nonlinear_structure_data_model |> mutate(type = "model"),
                       scaled_one_nonlinear_structure_data |> mutate(type = "data"))

## Set the maximum difference as the criteria
distance_df_small_edges_one_nonlinear_structure <- distance_one_nonlinear_structure |>
  dplyr::filter(distance < benchmark_one_nonlinear_structure)

# Visualize with langevitour
langevitour(df |> dplyr::select(-type),
            lineFrom = distance_df_small_edges_one_nonlinear_structure$from,
            lineTo = distance_df_small_edges_one_nonlinear_structure$to,
            group = df$type,
            pointSize = append(rep(1.5, NROW(scaled_one_nonlinear_structure_data_model)), rep(1, NROW(scaled_one_nonlinear_structure_data))),
            levelColors = c("#000000", "#33a02c"),
            lineColors = rep("#33a02c", nrow(distance_df_small_edges_one_nonlinear_structure)))


