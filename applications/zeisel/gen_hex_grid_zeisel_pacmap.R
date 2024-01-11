library(readr)
library(dplyr)

set.seed(20240110)

source("quollr_code.R", local = TRUE)
source("nldr_code.R", local = TRUE)

## Import data
df_zeisel <- read_rds("data/zeisel/zeisel.rds")
training_data_zeisel <- read_rds("data/zeisel/zeisel_training.rds")

PaCMAP_zeisel <- read_rds("data/zeisel/zeisel_pacmap.rds")

plot_PaCMAP_2D(PaCMAP_zeisel) + #ggtitle("(b)") +
  theme_linedraw() +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate(geom = 'text', label = 'b', x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 3)

num_bins_zeisel <- calculate_effective_x_bins(.data = PaCMAP_zeisel, x = PaCMAP1,
                                              cell_area = 1) ##31

shape_value_zeisel <- calculate_effective_shape_value(.data = PaCMAP_zeisel,
                                                      x = PaCMAP1, y = PaCMAP2) ## 0.8421752

## To extract bin centroids
hexbin_data_object_zeisel <- extract_hexbin_centroids(nldr_df = PaCMAP_zeisel,
                                                      num_bins = num_bins_zeisel,
                                                      shape_val = shape_value_zeisel, x = PaCMAP1, y = PaCMAP2)

df_bin_centroids_zeisel <- hexbin_data_object_zeisel$hexdf_data

##########

## Data set with all possible centroids in the hexagonal grid

full_centroid_df <- generate_full_grid_centroids(df_bin_centroids_zeisel)

## Generate all coordinates of hexagons
hex_grid <- full_hex_grid(full_centroid_df)

hex_full_count_df <- generate_full_grid_info(df_bin_centroids_zeisel)

write_rds(hex_full_count_df, "data/zeisel/zeisel_pacmap_hex_31.rds")

##########

ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_point(data = PaCMAP_zeisel, aes(x = PaCMAP1, y = PaCMAP2), alpha = 0.5) +
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


PaCMAP_zeisel_with_hb_id <- PaCMAP_zeisel |>
  dplyr::mutate(hb_id = hexbin_data_object_zeisel$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_zeisel <- dplyr::bind_cols(training_data_zeisel |> dplyr::select(-ID), PaCMAP_zeisel_with_hb_id)

## Averaged on high-D
df_bin_zeisel <- avg_highD_data(.data = df_all_zeisel, column_start_text = "PC")

## Triangulate bin centroids
tr1_object_zeisel <- triangulate_bin_centroids(df_bin_centroids_zeisel, x, y)
tr_from_to_df_zeisel <- generate_edge_info(triangular_object = tr1_object_zeisel)

## Compute 2D distances
distance_zeisel <- cal_2D_dist(.data = tr_from_to_df_zeisel)

## To find the benchmark value
benchmark_zeisel <- find_benchmark_value(.data = distance_zeisel, distance_col = distance)
##1.8

trimesh_zeisel_pacmap <- ggplot(df_bin_centroids_zeisel, aes(x = x, y = y)) +
  geom_segment(data = tr_from_to_df_zeisel, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
  geom_point(size = 2, colour = "#33a02c") +
  coord_equal()

# ggplot(df_bin_centroids, aes(x = x, y = y)) +
# geom_point(size = 1, colour = "#33a02c") +
# geom_trimesh() +
# coord_equal()

trimesh_zeisel_pacmap <- trimesh_zeisel_pacmap +
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

trimesh_gr_zeisel_pacmap <- colour_long_edges(.data = distance_zeisel, benchmark_value = benchmark_zeisel,
                                            triangular_object = tr1_object_zeisel, distance_col = distance)

trimesh_gr_zeisel_pacmap <- trimesh_gr_zeisel_pacmap +
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

trimesh_removed_zeisel_pacmap <- remove_long_edges(.data = distance_zeisel, benchmark_value = benchmark_zeisel,
                                                 triangular_object = tr1_object_zeisel, distance_col = distance)

trimesh_removed_zeisel_pacmap <- trimesh_removed_zeisel_pacmap +
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



tour1_zeisel_pacmap <- show_langevitour(df_all_zeisel, df_bin_zeisel,
                                      df_bin_centroids_zeisel, benchmark_value = benchmark_zeisel,
                                      distance = distance_zeisel, distance_col = distance, col_start = "PC")
