library(readr)
library(dplyr)

set.seed(20240110)
source("quollr_code.R", local = TRUE)

## Import data
data <- read_rds("data/s_curve/data_s_curve.rds")
training_data <- read_rds("data/s_curve/data_s_curve_training.rds")
test_data <- read_rds("data/s_curve/data_s_curve_test.rds")

tSNE_s_curve <- read_rds("data/s_curve/s_curve_tsne_27.rds")
UMAP_s_curve <- read_rds("data/s_curve/s_curve_umap.rds")
PHATE_s_curve <- read_rds("data/s_curve/s_curve_phate.rds")
TriMAP_s_curve <- read_rds("data/s_curve/s_curve_trimap.rds")
PaCMAP_s_curve <- read_rds("data/s_curve/s_curve_pacmap.rds")


## tSNE

num_bins_tsne_s_curve <- 8
shape_val_tsne_s_curve <- calculate_effective_shape_value(.data = tSNE_s_curve,
                                                          x = tSNE1, y = tSNE2) ## 1.259938
## To extract bin centroids
hexbin_data_object_tsne_s_curve <- extract_hexbin_centroids(nldr_df = tSNE_s_curve, num_bins = num_bins_tsne_s_curve, shape_val = shape_val_tsne_s_curve, x = tSNE1, y = tSNE2)

df_bin_centroids_tsne_s_curve <- hexbin_data_object_tsne_s_curve$hexdf_data

tSNE_data_with_hb_id_s_curve <- tSNE_s_curve |>
  dplyr::mutate(hb_id = hexbin_data_object_tsne_s_curve$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_tsne_s_curve <- dplyr::bind_cols(training_data |> dplyr::select(-ID), tSNE_data_with_hb_id_s_curve)

## Averaged on high-D
df_bin_tsne_s_curve <- avg_highD_data(.data = df_all_tsne_s_curve)

## Triangulate bin centroids
tr1_object_tsne_s_curve <- triangulate_bin_centroids(df_bin_centroids_tsne_s_curve, x, y)
tr_from_to_df_tsne_s_curve <- generate_edge_info(triangular_object = tr1_object_tsne_s_curve)

# ggplot(df_bin_centroids_tsne_s_curve, aes(x = x, y = y)) +
#   geom_segment(data = tr_from_to_df_tsne_s_curve, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
#   geom_point(size = 2, colour = "#33a02c") +
#   coord_equal()


## Compute 2D distances
distance_tsne_s_curve <- cal_2D_dist(.data = tr_from_to_df_tsne_s_curve)

## To find the benchmark value
benchmark_tsne_s_curve <- find_benchmark_value(.data = distance_tsne_s_curve, distance_col = distance)

# colour_long_edges(.data = distance_tsne_s_curve, benchmark_value = benchmark_tsne_s_curve,
#                   triangular_object = tr1_object_tsne_s_curve, distance_col = distance)

trimesh_removed_tsne_s_curve <- remove_long_edges(.data = distance_tsne_s_curve, benchmark_value = benchmark_tsne_s_curve,
                                                  triangular_object = tr1_object_tsne_s_curve, distance_col = distance)

trimesh_removed_tsne_s_curve <- trimesh_removed_tsne_s_curve +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))

tour_tsne_s_curve <- show_langevitour(df_all_tsne_s_curve, df_bin_tsne_s_curve, df_bin_centroids_tsne_s_curve, benchmark_value = benchmark_tsne_s_curve, distance = distance_tsne_s_curve, distance_col = distance)



## UMAP

num_bins_umap_s_curve <- 6
shape_val_umap_s_curve <- calculate_effective_shape_value(.data = UMAP_s_curve,
                                                          x = UMAP1, y = UMAP2) ## 1.259938
## To extract bin centroids
hexbin_data_object_umap_s_curve <- extract_hexbin_centroids(nldr_df = UMAP_s_curve, num_bins = num_bins_umap_s_curve, shape_val = shape_val_umap_s_curve, x = UMAP1, y = UMAP2)

df_bin_centroids_umap_s_curve <- hexbin_data_object_umap_s_curve$hexdf_data

UMAP_data_with_hb_id_s_curve <- UMAP_s_curve |>
  dplyr::mutate(hb_id = hexbin_data_object_umap_s_curve$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_umap_s_curve <- dplyr::bind_cols(training_data |> dplyr::select(-ID), UMAP_data_with_hb_id_s_curve)

## Averaged on high-D
df_bin_umap_s_curve <- avg_highD_data(.data = df_all_umap_s_curve)

## Triangulate bin centroids
tr1_object_umap_s_curve <- triangulate_bin_centroids(df_bin_centroids_umap_s_curve, x, y)
tr_from_to_df_umap_s_curve <- generate_edge_info(triangular_object = tr1_object_umap_s_curve)

# ggplot(df_bin_centroids_umap_s_curve, aes(x = x, y = y)) +
#   geom_segment(data = tr_from_to_df_umap_s_curve, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
#   geom_point(size = 2, colour = "#33a02c") +
#   coord_equal()


## Compute 2D distances
distance_umap_s_curve <- cal_2D_dist(.data = tr_from_to_df_umap_s_curve)

## To find the benchmark value
benchmark_umap_s_curve <- find_benchmark_value(.data = distance_umap_s_curve, distance_col = distance)

# colour_long_edges(.data = distance_umap_s_curve, benchmark_value = benchmark_umap_s_curve,
#                   triangular_object = tr1_object_umap_s_curve, distance_col = distance)

trimesh_removed_umap_s_curve <- remove_long_edges(.data = distance_umap_s_curve, benchmark_value = benchmark_umap_s_curve,
                                     triangular_object = tr1_object_umap_s_curve, distance_col = distance)

trimesh_removed_umap_s_curve <- trimesh_removed_umap_s_curve +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))

tour_umap_s_curve <- show_langevitour(df_all_umap_s_curve, df_bin_umap_s_curve, df_bin_centroids_umap_s_curve, benchmark_value = benchmark_umap_s_curve, distance = distance_umap_s_curve, distance_col = distance)


## TriMAP

num_bins_trimap_s_curve <- 6
shape_val_trimap_s_curve <- calculate_effective_shape_value(.data = TriMAP_s_curve,
                                                            x = TriMAP1, y = TriMAP2) ## 1.259938
## To extract bin centroids
hexbin_data_object_trimap_s_curve <- extract_hexbin_centroids(nldr_df = TriMAP_s_curve, num_bins = num_bins_trimap_s_curve, shape_val = shape_val_trimap_s_curve, x = TriMAP1, y = TriMAP2)

df_bin_centroids_trimap_s_curve <- hexbin_data_object_trimap_s_curve$hexdf_data

TriMAP_data_with_hb_id_s_curve <- TriMAP_s_curve |>
  dplyr::mutate(hb_id = hexbin_data_object_trimap_s_curve$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_trimap_s_curve <- dplyr::bind_cols(training_data |> dplyr::select(-ID), TriMAP_data_with_hb_id_s_curve)

## Averaged on high-D
df_bin_trimap_s_curve <- avg_highD_data(.data = df_all_trimap_s_curve)

## Triangulate bin centroids
tr1_object_trimap_s_curve <- triangulate_bin_centroids(df_bin_centroids_trimap_s_curve, x, y)
tr_from_to_df_trimap_s_curve <- generate_edge_info(triangular_object = tr1_object_trimap_s_curve)

# ggplot(df_bin_centroids_trimap_s_curve, aes(x = x, y = y)) +
#   geom_segment(data = tr_from_to_df_trimap_s_curve, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
#   geom_point(size = 2, colour = "#33a02c") +
#   coord_equal()


## Compute 2D distances
distance_trimap_s_curve <- cal_2D_dist(.data = tr_from_to_df_trimap_s_curve)

## To find the benchmark value
benchmark_trimap_s_curve <- find_benchmark_value(.data = distance_trimap_s_curve, distance_col = distance)

# colour_long_edges(.data = distance_trimap_s_curve, benchmark_value = benchmark_trimap_s_curve,
#                   triangular_object = tr1_object_trimap_s_curve, distance_col = distance)

trimesh_removed_trimap_s_curve <- remove_long_edges(.data = distance_trimap_s_curve, benchmark_value = benchmark_trimap_s_curve,
                                                    triangular_object = tr1_object_trimap_s_curve, distance_col = distance)

trimesh_removed_trimap_s_curve <- trimesh_removed_trimap_s_curve +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))

tour_trimap_s_curve <- show_langevitour(df_all_trimap_s_curve, df_bin_trimap_s_curve, df_bin_centroids_trimap_s_curve, benchmark_value = benchmark_trimap_s_curve, distance = distance_trimap_s_curve, distance_col = distance)



## PacMAP

num_bins_pacmap_s_curve <- 10
shape_val_pacmap_s_curve <- calculate_effective_shape_value(.data = PaCMAP_s_curve,
                                                            x = PaCMAP1, y = PaCMAP2) ## 1.259938
## To extract bin centroids
hexbin_data_object_pacmap_s_curve <- extract_hexbin_centroids(nldr_df = PaCMAP_s_curve, num_bins = num_bins_pacmap_s_curve, shape_val = shape_val_pacmap_s_curve, x = PaCMAP1, y = PaCMAP2)

df_bin_centroids_pacmap_s_curve <- hexbin_data_object_pacmap_s_curve$hexdf_data

PaCMAP_data_with_hb_id_s_curve <- PaCMAP_s_curve |>
  dplyr::mutate(hb_id = hexbin_data_object_pacmap_s_curve$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_pacmap_s_curve <- dplyr::bind_cols(training_data |> dplyr::select(-ID), PaCMAP_data_with_hb_id_s_curve)

## Averaged on high-D
df_bin_pacmap_s_curve <- avg_highD_data(.data = df_all_pacmap_s_curve)

## Triangulate bin centroids
tr1_object_pacmap_s_curve <- triangulate_bin_centroids(df_bin_centroids_pacmap_s_curve, x, y)
tr_from_to_df_pacmap_s_curve <- generate_edge_info(triangular_object = tr1_object_pacmap_s_curve)

# ggplot(df_bin_centroids_pacmap_s_curve, aes(x = x, y = y)) +
#   geom_segment(data = tr_from_to_df_pacmap_s_curve, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
#   geom_point(size = 2, colour = "#33a02c") +
#   coord_equal()


## Compute 2D distances
distance_pacmap_s_curve <- cal_2D_dist(.data = tr_from_to_df_pacmap_s_curve)

## To find the benchmark value
benchmark_pacmap_s_curve <- find_benchmark_value(.data = distance_pacmap_s_curve, distance_col = distance)

# colour_long_edges(.data = distance_pacmap_s_curve, benchmark_value = benchmark_pacmap_s_curve,
#                   triangular_object = tr1_object_pacmap_s_curve, distance_col = distance)

trimesh_removed_pacmap_s_curve <- remove_long_edges(.data = distance_pacmap_s_curve, benchmark_value = benchmark_pacmap_s_curve,
                                     triangular_object = tr1_object_pacmap_s_curve, distance_col = distance)

trimesh_removed_pacmap_s_curve <- trimesh_removed_pacmap_s_curve +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))

tour_pacmap_s_curve <- show_langevitour(df_all_pacmap_s_curve, df_bin_pacmap_s_curve, df_bin_centroids_pacmap_s_curve, benchmark_value = benchmark_pacmap_s_curve, distance = distance_pacmap_s_curve, distance_col = distance)

## PHATE

num_bins_phate_s_curve <- 8
shape_val_phate_s_curve <- calculate_effective_shape_value(.data = PHATE_s_curve,
                                                           x = PHATE1, y = PHATE2) ## 1.259938
## To extract bin centroids
hexbin_data_object_phate_s_curve <- extract_hexbin_centroids(nldr_df = PHATE_s_curve, num_bins = num_bins_phate_s_curve, shape_val = shape_val_phate_s_curve, x = PHATE1, y = PHATE2)

df_bin_centroids_phate_s_curve <- hexbin_data_object_phate_s_curve$hexdf_data

## Identify bins with low-density
identify_rm_bins <- find_low_density_hexagons(df_bin_centroids_phate_s_curve, num_bins_phate_s_curve, benchmark_rm_hex = 0.06)

df_bin_centroids_phate_s_curve <- df_bin_centroids_phate_s_curve |>
  filter(!(hexID %in% identify_rm_bins))

PHATE_data_with_hb_id_s_curve <- PHATE_s_curve |>
  dplyr::mutate(hb_id = hexbin_data_object_phate_s_curve$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_phate_s_curve <- dplyr::bind_cols(training_data |> dplyr::select(-ID), PHATE_data_with_hb_id_s_curve)

## Averaged on high-D
df_bin_phate_s_curve <- avg_highD_data(.data = df_all_phate_s_curve)

## Triangulate bin centroids
tr1_object_phate_s_curve <- triangulate_bin_centroids(df_bin_centroids_phate_s_curve, x, y)
tr_from_to_df_phate_s_curve <- generate_edge_info(triangular_object = tr1_object_phate_s_curve)

# ggplot(df_bin_centroids_phate_s_curve, aes(x = x, y = y)) +
#   geom_segment(data = tr_from_to_df_phate_s_curve, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
#   geom_point(size = 2, colour = "#33a02c") +
#   coord_equal()


## Compute 2D distances
distance_phate_s_curve <- cal_2D_dist(.data = tr_from_to_df_phate_s_curve)

## To find the benchmark value
benchmark_phate_s_curve <- find_benchmark_value(.data = distance_phate_s_curve, distance_col = distance)

# colour_long_edges(.data = distance_phate_s_curve, benchmark_value = benchmark_phate_s_curve,
#                   triangular_object = tr1_object_phate_s_curve, distance_col = distance)

trimesh_removed_phate_s_curve <- remove_long_edges(.data = distance_phate_s_curve, benchmark_value = benchmark_phate_s_curve,
                                                   triangular_object = tr1_object_phate_s_curve, distance_col = distance)

trimesh_removed_phate_s_curve <- trimesh_removed_phate_s_curve +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 7))

tour_phate_s_curve <- show_langevitour(df_all_phate_s_curve, df_bin_phate_s_curve, df_bin_centroids_phate_s_curve, benchmark_value = benchmark_phate_s_curve, distance = distance_phate_s_curve, distance_col = distance)
