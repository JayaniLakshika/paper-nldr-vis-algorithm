library(dplyr)
library(rsample)
library(readr)
library(umap)

set.seed(20240110)

source("quollr_code.R", local = TRUE)
source("nldr_code.R", local = TRUE)

#save(gene.list,label.list,umap.list,plots.list,file = "./results/pbmc3k_clustering_UMAP.RData")

data_pca <- read_rds(file = "data/pbmc/pbmc_3k_festem/pbmc_pca_50.rds")

data_pca <- data_pca[,1:15]

training_data_pbmc <- data_pca |>
  mutate(ID = row_number())

## UMAP

UMAP_pbmc <- read_rds("data/pbmc/pbmc_umap_5_min_dist_0.99_metric_cosine.rds")

#num_bins_pbmc <- calculate_effective_x_bins(.data = UMAP_pbmc, x = UMAP1,
#                                            cell_area = 1) ##23

shape_value_pbmc <- calculate_effective_shape_value(.data = UMAP_pbmc,
                                                    x = UMAP1, y = UMAP2) ## 0.8772751

# ## To extract bin centroids
# hexbin_data_object_pbmc <- extract_hexbin_centroids(nldr_df = UMAP_pbmc,
#                                                     num_bins = num_bins_pbmc,
#                                                     shape_val = shape_value_pbmc, x = UMAP1, y = UMAP2)
#
# df_bin_centroids_pbmc <- hexbin_data_object_pbmc$hexdf_data
#
# ##########
#
# ## Data set with all possible centroids in the hexagonal grid
#
# full_centroid_df <- generate_full_grid_centroids(df_bin_centroids_pbmc)
#
# ## Generate all coordinates of hexagons
# hex_grid <- full_hex_grid(full_centroid_df)
#
# hex_full_count_df <- generate_full_grid_info(df_bin_centroids_pbmc)

#write_rds(hex_full_count_df, "data/pbmc/pbmc_3k_festem/pbmc_umap_hex_23.rds")

##########
## Need only 83 non-empty bins
shape_val <- shape_value_pbmc
non_empty_bins <- 83

#non_empty_bins: number of non-empty bins needed

# find_non_empty_bins <- function(.data, x = "UMAP1", y = "UMAP2", shape_val, non_empty_bins) {
#
#   num_bins_x <- 1
#   ## To extract bin centroids
#   hexbin_data_object <- extract_hexbin_centroids(nldr_df = .data,
#                                                       num_bins = num_bins_x,
#                                                       shape_val = shape_val, x = x, y = y)
#   df_bin_centroids <- hexbin_data_object$hexdf_data
#
#   num_of_non_empty_bins <- df_bin_centroids$hexID |> length()
#
#   while (num_of_non_empty_bins < non_empty_bins) {
#
#     num_bins_x <- num_bins_x + 1
#
#
#     ## To extract bin centroids
#     hexbin_data_object <- extract_hexbin_centroids(nldr_df = .data,
#                                                         num_bins = num_bins_x,
#                                                         shape_val = shape_val, x = y, y = y)
#
#     df_bin_centroids <- hexbin_data_object$hexdf_data
#
#     num_of_non_empty_bins <- df_bin_centroids$hexID |> length()
#
#     if (num_of_non_empty_bins >= non_empty_bins) {
#
#       return(num_bins_x)
#       break
#
#     } else {
#       next
#
#     }
#
#   }
#
# }



num_bins_pbmc <- find_non_empty_bins(nldr_df = UMAP_pbmc, x = "UMAP1", y = "UMAP2", shape_val, non_empty_bins)


## To extract bin centroids
hexbin_data_object_pbmc <- extract_hexbin_centroids(nldr_df = UMAP_pbmc,
                                                    num_bins = num_bins_pbmc,
                                                    shape_val = shape_value_pbmc, x = UMAP1, y = UMAP2)

df_bin_centroids_pbmc <- hexbin_data_object_pbmc$hexdf_data

##########

## Data set with all possible centroids in the hexagonal grid

full_centroid_df <- generate_full_grid_centroids(df_bin_centroids_pbmc)

## Generate all coordinates of hexagons
hex_grid <- full_hex_grid(full_centroid_df)

hex_full_count_df <- generate_full_grid_info(df_bin_centroids_pbmc)

ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_text(aes(x = c_x, y = c_y, label = hexID)) +
  scale_fill_viridis_c(direction = -1, na.value = "#ffffff")

