library(readr)
library(dplyr)

set.seed(20240110)

source("quollr_code.R", local = TRUE)
source("nldr_code.R", local = TRUE)

## Import data
#df_kang <- read_rds("data/kang/kang.rds")
training_data_kang <- read_rds("data/kang/kang_pca_25.rds")
training_data_kang <- training_data_kang[, 1:25] |>
  mutate(ID = 1:NROW(training_data_kang))

tSNE_kang <- read_rds("data/kang/kang_tsne.rds")

num_bins_kang <- calculate_effective_x_bins(.data = tSNE_kang, x = tSNE1,
                                            cell_area = 1) ##78

num_bins_kang <- 30

shape_value_kang <- calculate_effective_shape_value(.data = tSNE_kang,
                                                    x = tSNE1, y = tSNE2) ## 0.8772751

## To extract bin centroids
hexbin_data_object_kang <- extract_hexbin_centroids(nldr_df = tSNE_kang,
                                                    num_bins = num_bins_kang,
                                                    shape_val = shape_value_kang, x = tSNE1, y = tSNE2)

df_bin_centroids_kang <- hexbin_data_object_kang$hexdf_data

##########

## Data set with all possible centroids in the hexagonal grid

full_centroid_df <- generate_full_grid_centroids(df_bin_centroids_kang)

## Generate all coordinates of hexagons
hex_grid <- full_hex_grid(full_centroid_df)

hex_full_count_df <- generate_full_grid_info(df_bin_centroids_kang)

write_rds(hex_full_count_df, "data/kang/kang_tsne_hex_78.rds")

##########

ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_point(data = tSNE_kang, aes(x = tSNE1, y = tSNE2), alpha = 0.5) +
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


tSNE_kang_with_hb_id <- tSNE_kang |>
  dplyr::mutate(hb_id = hexbin_data_object_kang$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_kang <- dplyr::bind_cols(training_data_kang |> dplyr::select(-ID), tSNE_kang_with_hb_id)

## Averaged on high-D
df_bin_kang <- avg_highD_data(.data = df_all_kang, column_start_text = "PC")

## Triangulate bin centroids
tr1_object_kang <- triangulate_bin_centroids(df_bin_centroids_kang, x, y)
tr_from_to_df_kang <- generate_edge_info(triangular_object = tr1_object_kang)

## Compute 2D distances
distance_kang <- cal_2D_dist(.data = tr_from_to_df_kang)

## To find the benchmark value
benchmark_kang <- find_benchmark_value(.data = distance_kang, distance_col = distance)
##1.8

trimesh_kang_tsne <- ggplot(df_bin_centroids_kang, aes(x = x, y = y)) +
  geom_segment(data = tr_from_to_df_kang, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
  geom_point(size = 2, colour = "#33a02c") +
  coord_equal()

# ggplot(df_bin_centroids, aes(x = x, y = y)) +
# geom_point(size = 1, colour = "#33a02c") +
# geom_trimesh() +
# coord_equal()

trimesh_kang_tsne <- trimesh_kang_tsne +
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

trimesh_gr_kang_tsne <- colour_long_edges(.data = distance_kang, benchmark_value = benchmark_kang,
                                          triangular_object = tr1_object_kang, distance_col = distance)

trimesh_gr_kang_tsne <- trimesh_gr_kang_tsne +
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

trimesh_removed_kang_tsne <- remove_long_edges(.data = distance_kang, benchmark_value = benchmark_kang,
                                               triangular_object = tr1_object_kang, distance_col = distance)

trimesh_removed_kang_tsne <- trimesh_removed_kang_tsne +
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



tour1_kang_tsne <- show_langevitour(df_all_kang, df_bin_kang,
                                    df_bin_centroids_kang, benchmark_value = benchmark_kang,
                                    distance = distance_kang, distance_col = distance, col_start = "PC")
