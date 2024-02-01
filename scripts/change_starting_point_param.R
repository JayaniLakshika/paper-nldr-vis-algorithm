library(readr)
# library(umap)
# library(dplyr)
# library(rsample)


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

ggplot(data = hex_grid, aes(x = x, y = y)) + geom_polygon(fill = "white", color = "black", aes(group = id)) +
  geom_point()

hex_full_count_df <- generate_full_grid_info(df_bin_centroids)

ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_text(aes(x = c_x, y = c_y, label = hexID)) +
  scale_fill_viridis_c(direction = -1, na.value = "#ffffff")

## Filter centroids with their hexIDs
hexbin_coord_all <- hex_full_count_df |>
  select(c_x, c_y, hexID) |>
  distinct()

hexbin_coord_all_new <- hexbin_coord_all |>
  mutate(c_x = c_x - (cell_diameter/2),
         c_y = c_y - (cell_diameter/2)) |>
  rename(c("x" = "c_x",
           "y" = "c_y")) |>
  select(-hexID)


## Generate all coordinates of hexagons
hex_grid_new <- full_hex_grid(hexbin_coord_all_new)

ggplot(data = hex_grid_new, aes(x = x, y = y)) + geom_polygon(fill = "white", color = "black", aes(group = id)) +
  geom_point()

ggplot(data = hex_grid_new, aes(x = x, y = y)) + geom_polygon(fill = "white", color = "black", aes(group = id)) +
  geom_point(data = hex_grid_new, aes(x = x, y = y), color = "black") +
  geom_point(data = hexbin_coord_all_new, aes(x = x, y = y), color = "black")

## Map hexID

vec1 <- stats::setNames(rep("", 2), c("x", "y"))  ## Define column names

full_grid_with_hexbin_id <- dplyr::bind_rows(vec1)[0, ]
full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
  dplyr::mutate_if(is.character, as.numeric)

for(i in 1:length(sort(unique(hex_grid_new$y)))) {

  ## Filter the data set with a specific y value
  specific_y_val_df <- hex_grid_new |>
    dplyr::filter(y == sort(unique(hex_grid_new$y))[i])

  ordered_x_df <- specific_y_val_df |>
    dplyr::arrange(x)

  full_grid_with_hexbin_id <- dplyr::bind_rows(full_grid_with_hexbin_id, ordered_x_df)
}

full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
  dplyr::mutate(hexID = row_number())

full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
  dplyr::rename("c_x" = "x",
                "c_y" = "y")

hex_full_count_df_new <- dplyr::bind_cols(hex_grid_new, full_grid_with_hexbin_id |> arrange(id) |> select(-id))

ggplot(data = hex_full_count_df_new, aes(x = x, y = y)) +
  geom_polygon(fill = "white", color = "black", aes(group = id)) +
  geom_text(data = full_grid_with_hexbin_id, aes(x = c_x, y = c_y, label = id))


#####



full_grid_with_polygon_id_df <- map_polygon_id(full_grid_with_hexbin_id, hex_grid_new)

full_grid_with_hexbin_id_rep <- full_grid_with_polygon_id_df |>
  dplyr::slice(rep(1:n(), each = 6)) |>
  dplyr::arrange(polygon_id)

hex_full_count_df_new <- dplyr::bind_cols(hex_grid_new, full_grid_with_hexbin_id_rep)




ggplot(data = hexbin_coord_all_new, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_text(aes(x = c_x, y = c_y, label = hexID)) +
  scale_fill_viridis_c(direction = -1, na.value = "#ffffff")

## To generate data set with point info
#full_grid_with_hexbin_id <- map_hexbin_id(full_centroid_df, df_bin_centroids)

## Add hexbin Id to 2D embeddings
UMAP_data_with_hb_id <- UMAP_data |>
  mutate(hb_id = hexbin_data_object$hb_data@cID)




#
# pts_df <- find_pts_in_hexbins(full_grid_with_hexbin_id, UMAP_data_with_hb_id)
#
# pts_df <- dplyr::full_join(full_grid_with_hexbin_id, pts_df, by = c("hexID" = "hexID")) |>
#   distinct()

## Find which point assign to which bin

## Select the point
UMAP_data_with_hb_id_spec <- UMAP_data_with_hb_id |>
  filter(row_number() == 1)

## Check it's within the same hexID bin or the nearest one
## If it's nearest one change the hexID to new one


