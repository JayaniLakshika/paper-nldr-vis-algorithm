library(readr)
library(dplyr)

set.seed(20240110)
source("quollr_code.R", local = TRUE)

## Import data
df_2 <- read_rds("data/five_gau_clusters/data_five_gau.rds")
training_data_1 <- read_rds("data/five_gau_clusters/data_five_gau_training.rds")
test_1 <- read_rds("data/five_gau_clusters/data_five_gau_test.rds")

tSNE_data_gau <- read_rds("data/five_gau_clusters/tsne_data_five_gau.rds")
UMAP_data_gau <- read_rds("data/five_gau_clusters/umap_data_five_gau.rds")
PHATE_data_gau <- read_rds("data/five_gau_clusters/phate_data_five_gau.rds")
TriMAP_data_gau <- read_rds("data/five_gau_clusters/trimap_data_five_gau.rds")
PaCMAP_data_gau <- read_rds("data/five_gau_clusters/pacmap_data_five_gau.rds")


### tSNE

num_bins_tsne_gau <- calculate_effective_x_bins(.data = tSNE_data_gau, x = tSNE1,
                                                cell_area = 1)
shape_val_tsne_gau <- calculate_effective_shape_value(.data = tSNE_data_gau,
                                                      x = tSNE1, y = tSNE2)
## To extract bin centroids
hexbin_data_object_tsne_gau <- extract_hexbin_centroids(nldr_df = tSNE_data_gau, num_bins = num_bins_tsne_gau, shape_val = shape_val_tsne_gau, x = tSNE1, y = tSNE2)

df_bin_centroids_tsne_gau <- hexbin_data_object_tsne_gau$hexdf_data

tSNE_data_with_hb_id_gau <- tSNE_data_gau |>
  dplyr::mutate(hb_id = hexbin_data_object_tsne_gau$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_tsne_gau <- dplyr::bind_cols(training_data_1 |> dplyr::select(-ID), tSNE_data_with_hb_id_gau)

## Averaged on high-D
df_bin_tsne_gau <- avg_highD_data(.data = df_all_tsne_gau)

## Triangulate bin centroids
tr1_object_tsne_gau <- triangulate_bin_centroids(df_bin_centroids_tsne_gau, x, y)
tr_from_to_df_tsne_gau <- generate_edge_info(triangular_object = tr1_object_tsne_gau)

## Compute 2D distances
distance_tsne_gau <- cal_2D_dist(.data = tr_from_to_df_tsne_gau)

## To find the benchmark value
benchmark_tsne_gau <- find_benchmark_value(.data = distance_tsne_gau, distance_col = distance)

tour_tsne_gau <- show_langevitour(df_all_tsne_gau, df_bin_tsne_gau, df_bin_centroids_tsne_gau, benchmark_value = benchmark_tsne_gau, distance = distance_tsne_gau, distance_col = distance)

### UMAP

num_bins_umap_gau <- calculate_effective_x_bins(.data = UMAP_data_gau, x = UMAP1,
                                                cell_area = 1)
shape_val_umap_gau <- calculate_effective_shape_value(.data = UMAP_data_gau,
                                                      x = UMAP1, y = UMAP2)
## To extract bin centroids
hexbin_data_object_umap_gau <- extract_hexbin_centroids(nldr_df = UMAP_data_gau, num_bins = num_bins_umap_gau, shape_val = shape_val_umap_gau, x = UMAP1, y = UMAP2)

df_bin_centroids_umap_gau <- hexbin_data_object_umap_gau$hexdf_data

UMAP_data_with_hb_id_gau <- UMAP_data_gau |>
  dplyr::mutate(hb_id = hexbin_data_object_umap_gau$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_umap_gau <- dplyr::bind_cols(training_data_1 |> dplyr::select(-ID), UMAP_data_with_hb_id_gau)

## Averaged on high-D
df_bin_umap_gau <- avg_highD_data(.data = df_all_umap_gau)

## Triangulate bin centroids
tr1_object_umap_gau <- triangulate_bin_centroids(df_bin_centroids_umap_gau, x, y)
tr_from_to_df_umap_gau <- generate_edge_info(triangular_object = tr1_object_umap_gau)

## Compute 2D distances
distance_umap_gau <- cal_2D_dist(.data = tr_from_to_df_umap_gau)

## To find the benchmark value
benchmark_umap_gau <- find_benchmark_value(.data = distance_umap_gau, distance_col = distance)

tour_umap_gau <- show_langevitour(df_all_umap_gau, df_bin_umap_gau, df_bin_centroids_umap_gau, benchmark_value = benchmark_umap_gau, distance = distance_umap_gau, distance_col = distance)


## PAHTE

num_bins_phate_gau <- calculate_effective_x_bins(.data = PHATE_data_gau, x = PHATE1,
                                                 cell_area = 1)
shape_val_phate_gau <- calculate_effective_shape_value(.data = PHATE_data_gau,
                                                       x = PHATE1, y = PHATE2)
## To extract bin centroids
hexbin_data_object_phate_gau <- extract_hexbin_centroids(nldr_df = PHATE_data_gau, num_bins = num_bins_phate_gau, shape_val = shape_val_phate_gau, x = PHATE1, y = PHATE2)

df_bin_centroids_phate_gau <- hexbin_data_object_phate_gau$hexdf_data

PHATE_data_with_hb_id_gau <- PHATE_data_gau |>
  dplyr::mutate(hb_id = hexbin_data_object_phate_gau$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_phate_gau <- dplyr::bind_cols(training_data_1 |> dplyr::select(-ID), PHATE_data_with_hb_id_gau)

## Averaged on high-D
df_bin_phate_gau <- avg_highD_data(.data = df_all_phate_gau)

## Triangulate bin centroids
tr1_object_phate_gau <- triangulate_bin_centroids(df_bin_centroids_phate_gau, x, y)
tr_from_to_df_phate_gau <- generate_edge_info(triangular_object = tr1_object_phate_gau)

## Compute 2D distances
distance_phate_gau <- cal_2D_dist(.data = tr_from_to_df_phate_gau)

## To find the benchmark value
benchmark_phate_gau <- find_benchmark_value(.data = distance_phate_gau, distance_col = distance)

tour_phate_gau <- show_langevitour(df_all_phate_gau, df_bin_phate_gau, df_bin_centroids_phate_gau, benchmark_value = benchmark_phate_gau, distance = distance_phate_gau, distance_col = distance)

## TriMAP

num_bins_trimap_gau <- calculate_effective_x_bins(.data = TriMAP_data_gau, x = TriMAP1,
                                                  cell_area = 1)
shape_val_trimap_gau <- calculate_effective_shape_value(.data = TriMAP_data_gau,
                                                        x = TriMAP1, y = TriMAP2)
## To extract bin centroids
hexbin_data_object_trimap_gau <- extract_hexbin_centroids(nldr_df = TriMAP_data_gau, num_bins = num_bins_trimap_gau, shape_val = shape_val_trimap_gau, x = TriMAP1, y = TriMAP2)

df_bin_centroids_trimap_gau <- hexbin_data_object_trimap_gau$hexdf_data

TriMAP_data_with_hb_id_gau <- TriMAP_data_gau |>
  dplyr::mutate(hb_id = hexbin_data_object_trimap_gau$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_trimap_gau <- dplyr::bind_cols(training_data_1 |> dplyr::select(-ID), TriMAP_data_with_hb_id_gau)

## Averaged on high-D
df_bin_trimap_gau <- avg_highD_data(.data = df_all_trimap_gau)

## Triangulate bin centroids
tr1_object_trimap_gau <- triangulate_bin_centroids(df_bin_centroids_trimap_gau, x, y)
tr_from_to_df_trimap_gau <- generate_edge_info(triangular_object = tr1_object_trimap_gau)

## Compute 2D distances
distance_trimap_gau <- cal_2D_dist(.data = tr_from_to_df_trimap_gau)

## To find the benchmark value
benchmark_trimap_gau <- find_benchmark_value(.data = distance_trimap_gau, distance_col = distance)

tour_trimap_gau <- show_langevitour(df_all_trimap_gau, df_bin_trimap_gau, df_bin_centroids_trimap_gau, benchmark_value = benchmark_trimap_gau, distance = distance_trimap_gau, distance_col = distance)


## PaCMAP

num_bins_pacmap_gau <- calculate_effective_x_bins(.data = PaCMAP_data_gau, x = PaCMAP1,
                                                  cell_area = 1)
shape_val_pacmap_gau <- calculate_effective_shape_value(.data = PaCMAP_data_gau,
                                                        x = PaCMAP1, y = PaCMAP2)
## To extract bin centroids
hexbin_data_object_pacmap_gau <- extract_hexbin_centroids(nldr_df = PaCMAP_data_gau, num_bins = num_bins_pacmap_gau, shape_val = shape_val_pacmap_gau, x = PaCMAP1, y = PaCMAP2)

df_bin_centroids_pacmap_gau <- hexbin_data_object_pacmap_gau$hexdf_data

PaCMAP_data_with_hb_id_gau <- PaCMAP_data_gau |>
  dplyr::mutate(hb_id = hexbin_data_object_pacmap_gau$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_pacmap_gau <- dplyr::bind_cols(training_data_1 |> dplyr::select(-ID), PaCMAP_data_with_hb_id_gau)

## Averaged on high-D
df_bin_pacmap_gau <- avg_highD_data(.data = df_all_pacmap_gau)

## Triangulate bin centroids
tr1_object_pacmap_gau <- triangulate_bin_centroids(df_bin_centroids_pacmap_gau, x, y)
tr_from_to_df_pacmap_gau <- generate_edge_info(triangular_object = tr1_object_pacmap_gau)

## Compute 2D distances
distance_pacmap_gau <- cal_2D_dist(.data = tr_from_to_df_pacmap_gau)

## To find the benchmark value
benchmark_pacmap_gau <- find_benchmark_value(.data = distance_pacmap_gau, distance_col = distance)

tour_pacmap_gau <- show_langevitour(df_all_pacmap_gau, df_bin_pacmap_gau, df_bin_centroids_pacmap_gau, benchmark_value = benchmark_pacmap_gau, distance = distance_pacmap_gau, distance_col = distance)
