library(readr)
# library(umap)
# library(dplyr)
# library(rsample)


set.seed(20240110)

source("quollr_code.R", local = TRUE)

## Read UMAP data
UMAP_data <- read_rds("data/s_curve/s_curve_umap.rds")


####num_bins_x = 3####################
num_bins <- 3
shape_val <- calculate_effective_shape_value(.data = UMAP_data, x = UMAP1, y = UMAP2)

hexbin_data_object_loop <- extract_hexbin_centroids(UMAP_data, num_bins, shape_val)

df_bin_centroids_loop <- hexbin_data_object_loop$hexdf_data

## Data set with all possible centroids in the hexagonal grid

full_centroid_df_loop <- generate_full_grid_centroids(df_bin_centroids_loop)

## To map hexID to hexbin centroids in the full grid

vec1 <- stats::setNames(rep("", 2), c("x", "y"))  ## Define column names

full_grid_with_hexbin_id_loop <- dplyr::bind_rows(vec1)[0, ]
full_grid_with_hexbin_id_loop <- full_grid_with_hexbin_id_loop |>
  dplyr::mutate_if(is.character, as.numeric)

for(i in 1:length(sort(unique(full_centroid_df_loop$y)))){

  ## Filter the data set with specific y value
  specific_y_val_df_loop <- full_centroid_df_loop |>
    dplyr::filter(y == sort(unique(full_centroid_df_loop$y))[i])

  ordered_x_df_loop <- specific_y_val_df_loop |>
    dplyr::arrange(x)

  full_grid_with_hexbin_id_loop <- dplyr::bind_rows(full_grid_with_hexbin_id_loop, ordered_x_df_loop)

}


full_grid_with_hexbin_id_loop <- full_grid_with_hexbin_id_loop |>
  dplyr::mutate(hexID = row_number())

full_grid_with_hexbin_id_loop <- full_grid_with_hexbin_id_loop |>
  dplyr::rename("c_x" = "x",
                "c_y" = "y")

full_grid_with_hexbin_id_loop <- dplyr::full_join(full_grid_with_hexbin_id_loop, df_bin_centroids_loop, by = c("hexID" = "hexID")) |>
  dplyr::select(-c(x, y))

full_grid_with_hexbin_id_loop <- full_grid_with_hexbin_id_loop |>
  dplyr::mutate(std_counts = counts/max(counts, na.rm = TRUE))

## Generate all coordinates of hexagons
hex_grid_loop <- full_hex_grid(full_centroid_df_loop)

full_grid_with_polygon_id_df_loop <- map_polygon_id(full_grid_with_hexbin_id_loop, hex_grid_loop)

full_grid_with_hexbin_id_rep_loop <- full_grid_with_polygon_id_df_loop |>
  dplyr::slice(rep(1:n(), each = 6)) |>
  dplyr::arrange(polygon_id)

hex_full_count_df_loop <- dplyr::bind_cols(hex_grid_loop, full_grid_with_hexbin_id_rep_loop)

## Run only once
write_rds(hex_full_count_df_loop, file = "data/s_curve/s_curve_hex_3.rds")

####num_bins_x = 8####################

num_bins <- 8
shape_val <- calculate_effective_shape_value(.data = UMAP_data, x = UMAP1, y = UMAP2)

hexbin_data_object_loop <- extract_hexbin_centroids(UMAP_data, num_bins, shape_val)

df_bin_centroids_loop <- hexbin_data_object_loop$hexdf_data

## Data set with all possible centroids in the hexagonal grid

full_centroid_df_loop <- generate_full_grid_centroids(df_bin_centroids_loop)

## To map hexID to hexbin centroids in the full grid

vec1 <- stats::setNames(rep("", 2), c("x", "y"))  ## Define column names

full_grid_with_hexbin_id_loop <- dplyr::bind_rows(vec1)[0, ]
full_grid_with_hexbin_id_loop <- full_grid_with_hexbin_id_loop |>
  dplyr::mutate_if(is.character, as.numeric)

for(i in 1:length(sort(unique(full_centroid_df_loop$y)))){

  ## Filter the data set with specific y value
  specific_y_val_df_loop <- full_centroid_df_loop |>
    dplyr::filter(y == sort(unique(full_centroid_df_loop$y))[i])

  ordered_x_df_loop <- specific_y_val_df_loop |>
    dplyr::arrange(x)

  full_grid_with_hexbin_id_loop <- dplyr::bind_rows(full_grid_with_hexbin_id_loop, ordered_x_df_loop)

}


full_grid_with_hexbin_id_loop <- full_grid_with_hexbin_id_loop |>
  dplyr::mutate(hexID = row_number())

full_grid_with_hexbin_id_loop <- full_grid_with_hexbin_id_loop |>
  dplyr::rename("c_x" = "x",
                "c_y" = "y")

full_grid_with_hexbin_id_loop <- dplyr::full_join(full_grid_with_hexbin_id_loop, df_bin_centroids_loop, by = c("hexID" = "hexID")) |>
  dplyr::select(-c(x, y))

full_grid_with_hexbin_id_loop <- full_grid_with_hexbin_id_loop |>
  dplyr::mutate(std_counts = counts/max(counts, na.rm = TRUE))

## Generate all coordinates of hexagons
hex_grid_loop <- full_hex_grid(full_centroid_df_loop)

full_grid_with_polygon_id_df_loop <- map_polygon_id(full_grid_with_hexbin_id_loop, hex_grid_loop)

full_grid_with_hexbin_id_rep_loop <- full_grid_with_polygon_id_df_loop |>
  dplyr::slice(rep(1:n(), each = 6)) |>
  dplyr::arrange(polygon_id)

hex_full_count_df_loop <- dplyr::bind_cols(hex_grid_loop, full_grid_with_hexbin_id_rep_loop)

## Run only once
write_rds(hex_full_count_df_loop, file = "data/s_curve/s_curve_hex_8.rds")


####num_bins_x = 15####################

num_bins <- 15
shape_val <- calculate_effective_shape_value(.data = UMAP_data, x = UMAP1, y = UMAP2)

hexbin_data_object_loop <- extract_hexbin_centroids(UMAP_data, num_bins, shape_val)

df_bin_centroids_loop <- hexbin_data_object_loop$hexdf_data

## Data set with all possible centroids in the hexagonal grid

full_centroid_df_loop <- generate_full_grid_centroids(df_bin_centroids_loop)

## To map hexID to hexbin centroids in the full grid

vec1 <- stats::setNames(rep("", 2), c("x", "y"))  ## Define column names

full_grid_with_hexbin_id_loop <- dplyr::bind_rows(vec1)[0, ]
full_grid_with_hexbin_id_loop <- full_grid_with_hexbin_id_loop |>
  dplyr::mutate_if(is.character, as.numeric)

for(i in 1:length(sort(unique(full_centroid_df_loop$y)))){

  ## Filter the data set with specific y value
  specific_y_val_df_loop <- full_centroid_df_loop |>
    dplyr::filter(y == sort(unique(full_centroid_df_loop$y))[i])

  ordered_x_df_loop <- specific_y_val_df_loop |>
    dplyr::arrange(x)

  full_grid_with_hexbin_id_loop <- dplyr::bind_rows(full_grid_with_hexbin_id_loop, ordered_x_df_loop)

}


full_grid_with_hexbin_id_loop <- full_grid_with_hexbin_id_loop |>
  dplyr::mutate(hexID = row_number())

full_grid_with_hexbin_id_loop <- full_grid_with_hexbin_id_loop |>
  dplyr::rename("c_x" = "x",
                "c_y" = "y")

full_grid_with_hexbin_id_loop <- dplyr::full_join(full_grid_with_hexbin_id_loop, df_bin_centroids_loop, by = c("hexID" = "hexID")) |>
  dplyr::select(-c(x, y))

full_grid_with_hexbin_id_loop <- full_grid_with_hexbin_id_loop |>
  dplyr::mutate(std_counts = counts/max(counts, na.rm = TRUE))

## Generate all coordinates of hexagons
hex_grid_loop <- full_hex_grid(full_centroid_df_loop)

full_grid_with_polygon_id_df_loop <- map_polygon_id(full_grid_with_hexbin_id_loop, hex_grid_loop)

full_grid_with_hexbin_id_rep_loop <- full_grid_with_polygon_id_df_loop |>
  dplyr::slice(rep(1:n(), each = 6)) |>
  dplyr::arrange(polygon_id)

hex_full_count_df_loop <- dplyr::bind_cols(hex_grid_loop, full_grid_with_hexbin_id_rep_loop)

## Run only once
write_rds(hex_full_count_df_loop, file = "data/s_curve/s_curve_hex_15.rds")

