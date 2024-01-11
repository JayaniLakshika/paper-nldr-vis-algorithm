library(readr)
library(dplyr)
library(ggplot2)

set.seed(20240110)

source("quollr_code.R", local = TRUE)
source("nldr_code.R", local = TRUE)

## Import data
df_s_curve <- read_csv("data/s_curve.csv")
training_data_s_curve <- read_rds("data/s_curve/s_curve_training.rds")

UMAP_s_curve <- read_rds("data/s_curve/s_curve_umap_7.rds")

plot_UMAP_2D(UMAP_s_curve) + #ggtitle("(b)") +
  theme_linedraw() +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate(geom = 'text', label = 'b', x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 3)

num_bins_s_curve <- calculate_effective_x_bins(.data = UMAP_s_curve, x = UMAP1,
                                              cell_area = 1) ##18

shape_value_s_curve <- calculate_effective_shape_value(.data = UMAP_s_curve,
                                                      x = UMAP1, y = UMAP2) ## 1.259938

## To extract bin centroids
hexbin_data_object_s_curve <- extract_hexbin_centroids(nldr_df = UMAP_s_curve,
                                                      num_bins = num_bins_s_curve,
                                                      shape_val = shape_value_s_curve, x = UMAP1, y = UMAP2)

df_bin_centroids_s_curve <- hexbin_data_object_s_curve$hexdf_data

##########

## Data set with all possible centroids in the hexagonal grid

full_centroid_df <- generate_full_grid_centroids(df_bin_centroids_s_curve)

## Generate all coordinates of hexagons
hex_grid <- full_hex_grid(full_centroid_df)

hex_full_count_df <- generate_full_grid_info(df_bin_centroids_s_curve)

write_rds(hex_full_count_df, "data/s_curve/s_curve_umap_hex_18.rds")

##########

ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_point(data = UMAP_s_curve, aes(x = UMAP1, y = UMAP2), alpha = 0.5) +
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


UMAP_s_curve_with_hb_id <- UMAP_s_curve |>
  dplyr::mutate(hb_id = hexbin_data_object_s_curve$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_s_curve <- dplyr::bind_cols(training_data_s_curve |> dplyr::select(-ID), UMAP_s_curve_with_hb_id)

## Averaged on high-D
df_bin_s_curve <- avg_highD_data(.data = df_all_s_curve, column_start_text = "PC")

## Triangulate bin centroids
tr1_object_s_curve <- triangulate_bin_centroids(df_bin_centroids_s_curve, x, y)
tr_from_to_df_s_curve <- generate_edge_info(triangular_object = tr1_object_s_curve)

## Compute 2D distances
distance_s_curve <- cal_2D_dist(.data = tr_from_to_df_s_curve)

## To find the benchmark value
benchmark_s_curve <- find_benchmark_value(.data = distance_s_curve, distance_col = distance)
##1.8

trimesh_s_curve_umap <- ggplot(df_bin_centroids_s_curve, aes(x = x, y = y)) +
  geom_segment(data = tr_from_to_df_s_curve, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
  geom_point(size = 2, colour = "#33a02c") +
  coord_equal()

# ggplot(df_bin_centroids, aes(x = x, y = y)) +
# geom_point(size = 1, colour = "#33a02c") +
# geom_trimesh() +
# coord_equal()

trimesh_s_curve_umap <- trimesh_s_curve_umap +
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

trimesh_gr_s_curve_umap <- colour_long_edges(.data = distance_s_curve, benchmark_value = benchmark_s_curve,
                                            triangular_object = tr1_object_s_curve, distance_col = distance)

trimesh_gr_s_curve_umap <- trimesh_gr_s_curve_umap +
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

trimesh_removed_s_curve_umap <- remove_long_edges(.data = distance_s_curve, benchmark_value = benchmark_s_curve,
                                                 triangular_object = tr1_object_s_curve, distance_col = distance)

trimesh_removed_s_curve_umap <- trimesh_removed_s_curve_umap +
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



tour1_s_curve_umap <- show_langevitour(df_all_s_curve, df_bin_s_curve,
                                      df_bin_centroids_s_curve, benchmark_value = benchmark_s_curve,
                                      distance = distance_s_curve, distance_col = distance, col_start = "PC")
