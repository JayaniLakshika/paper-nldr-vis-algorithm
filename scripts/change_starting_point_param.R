library(readr)
library(ggplot2)


set.seed(20240110)

source("quollr_code.R", local = TRUE)

## Read UMAP data
UMAP_data <- read_rds("data/s_curve/s_curve_umap.rds")


####num_bins_x = 8####################
num_bins <- 8

cell_diameter <- sqrt(2 * 1 / sqrt(3))
cell_diameter

cell_diameter/2

shape_val <- calculate_effective_shape_value(.data = UMAP_data, x = UMAP1, y = UMAP2)

hexbin_data_object <- extract_hexbin_centroids(UMAP_data, num_bins, shape_val)

df_bin_centroids <- hexbin_data_object$hexdf_data

## Data set with all possible centroids in the hexagonal grid

full_centroid_df <- generate_full_grid_centroids(df_bin_centroids)

## Generate all coordinates of hexagons
hex_grid <- full_hex_grid(full_centroid_df)

ggplot(data = hex_grid, aes(x = x, y = y)) + geom_polygon(fill = "white",
                                                          color = "black", aes(group = id)) +
  geom_point()

hex_full_count_df <- generate_full_grid_info(df_bin_centroids)
#write_rds(hex_full_count_df, file = "data/s_curve/s_curve_hex_8.rds")

ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_text(aes(x = c_x, y = c_y, label = hexID)) +
  scale_fill_viridis_c(direction = -1, na.value = "#ffffff")

## Filter centroids with their hexIDs
hexbin_coord_all <- hex_full_count_df |>
  dplyr::select(c_x, c_y, hexID) |>
  dplyr::distinct()

hexbin_coord_all_new <- hexbin_coord_all |>
  dplyr::mutate(c_x = c_x - (cell_diameter/2),
         c_y = c_y - (cell_diameter/2)) |>
  dplyr::rename(c("x" = "c_x",
           "y" = "c_y"))

ggplot(data = hex_full_count_df, aes(x = c_x, y = c_y)) +
  geom_point(color = "black") +
  geom_point(data = hexbin_coord_all_new, aes(x = x, y = y), color = "red")

## Generate all coordinates of hexagons
hex_grid_new <- full_hex_grid(hexbin_coord_all_new)

ggplot(data = hex_grid_new, aes(x = x, y = y)) + geom_polygon(fill = "white", color = "black", aes(group = id)) +
  geom_point()

ggplot(data = hex_grid_new, aes(x = x, y = y)) + geom_polygon(fill = "white", color = "black", aes(group = id)) +
  #geom_point(data = hex_grid_new, aes(x = x, y = y), color = "black") +
  geom_point(data = hexbin_coord_all_new, aes(x = x, y = y), color = "red") +
  geom_point(data = hex_full_count_df, aes(x = c_x, y = c_y), color = "blue")


# ## Map hexID
#
# vec1 <- stats::setNames(rep("", 2), c("x", "y"))  ## Define column names
#
# full_grid_with_hexbin_id <- dplyr::bind_rows(vec1)[0, ]
# full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
#   dplyr::mutate_if(is.character, as.numeric)
#
# for(i in 1:length(sort(unique(hex_grid_new$y)))) {
#
#   ## Filter the data set with a specific y value
#   specific_y_val_df <- hex_grid_new |>
#     dplyr::filter(y == sort(unique(hex_grid_new$y))[i])
#
#   ordered_x_df <- specific_y_val_df |>
#     dplyr::arrange(x)
#
#   full_grid_with_hexbin_id <- dplyr::bind_rows(full_grid_with_hexbin_id, ordered_x_df)
# }
#
# full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
#   dplyr::mutate(hexID = row_number())
#
# full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
#   dplyr::rename("c_x" = "x",
#                 "c_y" = "y")

hexbin_coord_all_new <- hexbin_coord_all_new |>
  dplyr::rename(c("c_x" = "x",
           "c_y" = "y"))

## Map the polygon ID to the hexagon coordinates
full_grid_with_polygon_id_df <- map_polygon_id(hexbin_coord_all_new, hex_grid_new)

full_grid_with_hexbin_id_rep <- full_grid_with_polygon_id_df |>
  dplyr::slice(rep(1:dplyr::n(), each = 6)) |>
  dplyr::arrange(polygon_id)

## Generate the dataset with polygon, and hexagon bin centroid coordinates
hex_full_count_df_new <- dplyr::bind_cols(hex_grid_new, full_grid_with_hexbin_id_rep)

ggplot(data = hex_full_count_df_new, aes(x = x, y = y)) +
  geom_polygon(fill = "white", color = "black", aes(group = polygon_id)) +
  geom_text(aes(x = c_x, y = c_y, label = hexID))

ggplot(data = hex_full_count_df_new, aes(x = x, y = y)) +
  geom_polygon(fill = "white", color = "black", aes(group = polygon_id)) +
  geom_polygon(data = hex_full_count_df, aes(x = x, y = y, group = polygon_id),
               fill = "white", color = "#feb24c")

##### Find counts within hexagons

## Add hexbin Id to 2D embeddings
UMAP_data_with_hb_id <- UMAP_data |>
  dplyr::mutate(hb_id = hexbin_data_object$hb_data@cID)

## Find which point assign to which bin

num_bins_x <- 8

nldr_df_with_new_hexID <- data.frame(matrix(ncol = 0, nrow = 0))

for (i in 1:NROW(UMAP_data_with_hb_id)) {

  ## Select the point
  UMAP_data_with_hb_id_spec <- UMAP_data_with_hb_id |>
    dplyr::filter(dplyr::row_number() == i)

  df_bin_centroids_coordinates_spec_bin_near1 <- hexbin_coord_all_new |>
    dplyr::filter((hexID == UMAP_data_with_hb_id_spec$hb_id[1]) |(hexID == (UMAP_data_with_hb_id_spec$hb_id[1] + (num_bins_x + 1))) | (hexID == (UMAP_data_with_hb_id_spec$hb_id[1] + num_bins_x)) | (hexID == (UMAP_data_with_hb_id_spec$hb_id[1] - (num_bins_x + 1))) | (hexID == (UMAP_data_with_hb_id_spec$hb_id[1] - num_bins_x)))

  UMAP_data_with_hb_id_spec <- UMAP_data_with_hb_id_spec |>
    dplyr::select(-ID) |>
    dplyr::rename("x" = "UMAP1",
           "y" = "UMAP2")

  df_bin_centroids_coordinates_spec_bin_near1 <- df_bin_centroids_coordinates_spec_bin_near1 |>
    dplyr::rename("x" = "c_x",
           "y" = "c_y",
           "hb_id" = "hexID")

  near_df_1 <- dplyr::bind_rows(UMAP_data_with_hb_id_spec, df_bin_centroids_coordinates_spec_bin_near1)

  near_df_1$distance <- lapply(seq(nrow(near_df_1)), function(x) {
    start <- unlist(near_df_1[1, c("x","y")])
    end <- unlist(near_df_1[x, c("x","y")])
    sqrt(sum((start - end)^2))})

  near_df_1$distance <- unlist(near_df_1$distance)

  near_df_1 <- near_df_1 |>
    dplyr::filter(dplyr::row_number() != 1) |>
    dplyr::arrange(distance)

  UMAP_data_with_hb_id_spec <- UMAP_data_with_hb_id_spec |>
    dplyr::select(-hb_id) |>
    dplyr::mutate(hb_id = near_df_1$hb_id[1])

  nldr_df_with_new_hexID <- dplyr::bind_rows(nldr_df_with_new_hexID, UMAP_data_with_hb_id_spec)

}


## Find counts within each hexagon

hb_id_with_counts <- nldr_df_with_new_hexID |>
  dplyr::count(hb_id) |>
  dplyr::mutate(counts = n,
                std_counts = n/max(n)) |>
  dplyr::select(-n)

hex_full_count_df_new <- dplyr::left_join(hex_full_count_df_new, hb_id_with_counts,
                                          by = c("hexID" = "hb_id"))

#write_rds(hex_full_count_df_new, file = "data/s_curve/s_curve_hex_8_shifted_hex_coord_df.rds")

ggplot(data = hex_full_count_df_new, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_text(aes(x = c_x, y = c_y, label = hexID)) +
  scale_fill_viridis_c(direction = -1, na.value = "#ffffff")

