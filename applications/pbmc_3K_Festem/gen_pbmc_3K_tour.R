library(readr)
library(dplyr)

set.seed(20240110)

source("quollr_code.R", local = TRUE)
source("nldr_code.R", local = TRUE)

## Import data
#df_pbmc <- read_rds("data/pbmc/pbmc.rds")
training_data_pbmc <- read_rds("data/pbmc/pbmc_3k_festem/pbmc_pca_50.rds")
training_data_pbmc <- training_data_pbmc[, 1:15] |>
  mutate(ID = 1:NROW(training_data_pbmc))

UMAP_pbmc <- read_rds("data/pbmc/pbmc_3k_festem/pbmc_umap.rds")

num_bins_pbmc <- calculate_effective_x_bins(.data = UMAP_pbmc, x = UMAP1,
                                              cell_area = 1) ##23

shape_value_pbmc <- calculate_effective_shape_value(.data = UMAP_pbmc,
                                                      x = UMAP1, y = UMAP2) ## 0.8772751

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

write_rds(hex_full_count_df, "data/pbmc/pbmc_3k_festem/pbmc_umap_hex_23.rds")

##########

ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_point(data = UMAP_pbmc, aes(x = UMAP1, y = UMAP2), alpha = 0.5) +
  scale_fill_viridis_c(direction = -1, na.value = "#ffffff", option = "C") +
  xlim(-5, 7) +
  ylim(-10, 10) +
  theme_void() +
  theme(legend.position="none", legend.direction="horizontal", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6)) +
  guides(fill = guide_colourbar(title = "Standardized count")) +
  annotate(geom = 'text', label = "a", x = -Inf, y = Inf, hjust = -0.3, vjust = 1, size = 3)


UMAP_pbmc_with_hb_id <- UMAP_pbmc |>
  dplyr::mutate(hb_id = hexbin_data_object_pbmc$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_pbmc <- dplyr::bind_cols(training_data_pbmc |> dplyr::select(-ID), UMAP_pbmc_with_hb_id)

## Averaged on high-D
df_bin_pbmc <- avg_highD_data(.data = df_all_pbmc, column_start_text = "PC")

## Triangulate bin centroids
tr1_object_pbmc <- triangulate_bin_centroids(df_bin_centroids_pbmc, x, y)
tr_from_to_df_pbmc <- generate_edge_info(triangular_object = tr1_object_pbmc)

## Compute 2D distances
distance_pbmc <- cal_2D_dist(.data = tr_from_to_df_pbmc)

## To find the benchmark value
benchmark_pbmc <- find_benchmark_value(.data = distance_pbmc, distance_col = distance)
##1.8

trimesh_pbmc_umap <- ggplot(df_bin_centroids_pbmc, aes(x = x, y = y)) +
  geom_segment(data = tr_from_to_df_pbmc, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
  geom_point(size = 2, colour = "#33a02c") +
  coord_equal()

# ggplot(df_bin_centroids, aes(x = x, y = y)) +
# geom_point(size = 1, colour = "#33a02c") +
# geom_trimesh() +
# coord_equal()

trimesh_pbmc_umap <- trimesh_pbmc_umap +
  #ggtitle("(a)") +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme_light() +
  theme(legend.position = "none", plot.title = element_text(size = 5, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()#change legend key width
  )
# theme(axis.text = element_text(size = 5),
#       axis.title = element_text(size = 7))

trimesh_gr_pbmc_umap <- colour_long_edges(.data = distance_pbmc, benchmark_value = benchmark_pbmc,
                                            triangular_object = tr1_object_pbmc, distance_col = distance)

trimesh_gr_pbmc_umap <- trimesh_gr_pbmc_umap +
  geom_point(size = 2, colour = "#33a02c") +
  #ggtitle("(b)") +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme_light() +
  #coord_equal() +
  theme(legend.position = "none", plot.title = element_text(size = 5, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()#change legend key width
  )

trimesh_removed_pbmc_umap <- remove_long_edges(.data = distance_pbmc, benchmark_value = benchmark_pbmc,
                                                 triangular_object = tr1_object_pbmc, distance_col = distance)

trimesh_removed_pbmc_umap <- trimesh_removed_pbmc_umap +
  geom_point(size = 2, colour = "#33a02c") +
  #ggtitle("(b)") +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme_light() +
  #coord_equal() +
  theme(legend.position = "none", plot.title = element_text(size = 5, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()#change legend key width
  )



tour1_pbmc_umap <- show_langevitour(df_all_pbmc, df_bin_pbmc,
                                      df_bin_centroids_pbmc, benchmark_value = benchmark_pbmc,
                                      distance = distance_pbmc, distance_col = distance, col_start = "PC")
