library(dplyr)
library(snedata)
library(ggflowchart)
library(purrr) ## map function
library(gridExtra) ## for grid.arrange
library(rsample)
library(DT)
library(ggbeeswarm)
library(ggplot2)
library(readr)

library(Rtsne)
library(umap)
library(phateR)
library(reticulate)
library(patchwork)

library(grid)
set.seed(20230531)

source(paste0(here::here(), "/paper-nldr-vis-algorithm/quollr_code.R"))
source(paste0(here::here(), "/paper-nldr-vis-algorithm/nldr_code.R"))

data <- read_csv(paste0(here::here(), "/paper-nldr-vis-algorithm/data/s_curve.csv"))

data <- data |>
  mutate(ID = row_number())

data_split <- initial_split(data)
training_data <- training(data_split) |>
  arrange(ID)
test_data <- testing(data_split) |>
  arrange(ID)

UMAP_fit <- umap(training_data |> dplyr::select(-ID), n_neighbors = 50, n_components =  2)

UMAP_data <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data)[1:(ncol(UMAP_data))] <- paste0(rep("UMAP",(ncol(UMAP_data))), 1:(ncol(UMAP_data)))

UMAP_data <- UMAP_data |>
  mutate(ID = training_data$ID)

num_bins_x <- calculate_effective_x_bins(.data = UMAP_data, x = UMAP1,
                                         cell_area = 1)

shape_value <- calculate_effective_shape_value(.data = UMAP_data,
                                               x = UMAP1, y = UMAP2)

## To extract bin centroids
hexbin_data_object <- extract_hexbin_centroids(nldr_df = UMAP_data, num_bins = num_bins_x, shape_val = shape_value)

df_bin_centroids <- hexbin_data_object$hexdf_data

df_bin_centroids_all <- df_bin_centroids

## Benchmark value to remove low-density hexagons is 0.255
identify_rm_bins <- find_low_density_hexagons(df_bin_centroids, num_bins_x)

df_bin_centroids <- df_bin_centroids |>
  filter(!(hexID %in% identify_rm_bins))

UMAP_data_with_hb_id <- UMAP_data |>
  dplyr::mutate(hb_id = hexbin_data_object$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all <- dplyr::bind_cols(training_data |> dplyr::select(-ID), UMAP_data_with_hb_id)

## Averaged on high-D
df_bin <- avg_highD_data(.data = df_all)

## Triangulate bin centroids
tr1_object <- triangulate_bin_centroids(df_bin_centroids, x, y)
tr_from_to_df <- generate_edge_info(triangular_object = tr1_object)

## Compute 2D distances
distance <- cal_2D_dist(.data = tr_from_to_df)

## To find the benchmark value
benchmark <- find_benchmark_value(.data = distance, distance_col = distance)

## With benchmark value 1.9
trimesh <- ggplot(df_bin_centroids, aes(x = x, y = y)) +
  geom_segment(data = tr_from_to_df, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
  geom_point(size = 2, colour = "#33a02c") +
  coord_equal()

trimesh <- trimesh +
  #ggtitle("(a)") +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme_light() +
  theme(legend.position = "none", plot.title = element_text(size = 5, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()#change legend key width
  )

trimesh_gr <- colour_long_edges(.data = distance, benchmark_value = benchmark,
                                triangular_object = tr1_object, distance_col = distance)

trimesh_gr <- trimesh_gr +
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

trimesh_removed <- remove_long_edges(.data = distance, benchmark_value = benchmark,
                                     triangular_object = tr1_object, distance_col = distance)

trimesh_removed <- trimesh_removed +
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



tour1 <- show_langevitour(df_all, df_bin, df_bin_centroids, benchmark_value = benchmark, distance = distance, distance_col = distance)

## To plot the distribution of distance
plot_dist <- function(distance_df){
  distance_df$group <- "1"
  dist_plot <- ggplot(distance_df, aes(x = group, y = distance)) +
    geom_quasirandom()+
    ylim(0, max(unlist(distance_df$distance))+ 0.5) + coord_flip()
  return(dist_plot)
}

distance_plot <- plot_dist(distance) +
  #ggtitle("(b)" ) +
  ylab(expression(d^{(2)})) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 12))

dist_p <- plot_dist(distance)  +
  geom_hline(yintercept = benchmark, colour = "red") +
  ylab(expression(d^{(2)})) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 12))

dist_p

trimesh + trimesh_gr +
  plot_layout(ncol = 2)



## Benchmark value to remove low-density hexagons is 0.1
df_bin_centroids <- df_bin_centroids_all

identify_rm_bins <- find_low_density_hexagons(df_bin_centroids, num_bins_x, benchmark_rm_hex = 0.2)

df_bin_centroids <- df_bin_centroids |>
  filter(!(hexID %in% identify_rm_bins))

UMAP_data_with_hb_id <- UMAP_data |>
  dplyr::mutate(hb_id = hexbin_data_object$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all <- dplyr::bind_cols(training_data |> dplyr::select(-ID), UMAP_data_with_hb_id)

## Averaged on high-D
df_bin <- avg_highD_data(.data = df_all)

## Triangulate bin centroids
tr1_object <- triangulate_bin_centroids(df_bin_centroids, x, y)
tr_from_to_df <- generate_edge_info(triangular_object = tr1_object)

## Compute 2D distances
distance <- cal_2D_dist(.data = tr_from_to_df)

## To find the benchmark value
benchmark <- find_benchmark_value(.data = distance, distance_col = distance)

## With benchmark value 1.9
trimesh <- ggplot(df_bin_centroids, aes(x = x, y = y)) +
  geom_segment(data = tr_from_to_df, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
  geom_point(size = 2, colour = "#33a02c") +
  coord_equal()

trimesh <- trimesh +
  #ggtitle("(a)") +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme_light() +
  theme(legend.position = "none", plot.title = element_text(size = 5, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()#change legend key width
  )

trimesh_gr <- colour_long_edges(.data = distance, benchmark_value = benchmark,
                                triangular_object = tr1_object, distance_col = distance)

trimesh_gr <- trimesh_gr +
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

trimesh_removed <- remove_long_edges(.data = distance, benchmark_value = benchmark,
                                     triangular_object = tr1_object, distance_col = distance)

trimesh_removed <- trimesh_removed +
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



tour1 <- show_langevitour(df_all, df_bin, df_bin_centroids, benchmark_value = benchmark, distance = distance, distance_col = distance)

## To plot the distribution of distance
plot_dist <- function(distance_df){
  distance_df$group <- "1"
  dist_plot <- ggplot(distance_df, aes(x = group, y = distance)) +
    geom_quasirandom()+
    ylim(0, max(unlist(distance_df$distance))+ 0.5) + coord_flip()
  return(dist_plot)
}

distance_plot <- plot_dist(distance) +
  #ggtitle("(b)" ) +
  ylab(expression(d^{(2)})) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 12))

dist_p <- plot_dist(distance)  +
  geom_hline(yintercept = benchmark, colour = "red") +
  ylab(expression(d^{(2)})) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 12))

dist_p

trimesh + trimesh_gr +
  plot_layout(ncol = 2)

## Benchmark value to remove low-density hexagons is 0.3
df_bin_centroids <- df_bin_centroids_all

identify_rm_bins <- find_low_density_hexagons(df_bin_centroids, num_bins_x, benchmark_rm_hex = 0.3)

df_bin_centroids <- df_bin_centroids |>
  filter(!(hexID %in% identify_rm_bins))

UMAP_data_with_hb_id <- UMAP_data |>
  dplyr::mutate(hb_id = hexbin_data_object$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all <- dplyr::bind_cols(training_data |> dplyr::select(-ID), UMAP_data_with_hb_id)

## Averaged on high-D
df_bin <- avg_highD_data(.data = df_all)

## Triangulate bin centroids
tr1_object <- triangulate_bin_centroids(df_bin_centroids, x, y)
tr_from_to_df <- generate_edge_info(triangular_object = tr1_object)

## Compute 2D distances
distance <- cal_2D_dist(.data = tr_from_to_df)

## To find the benchmark value
benchmark <- find_benchmark_value(.data = distance, distance_col = distance)

## With benchmark value 1.9
trimesh <- ggplot(df_bin_centroids, aes(x = x, y = y)) +
  geom_segment(data = tr_from_to_df, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
  geom_point(size = 2, colour = "#33a02c") +
  coord_equal()

trimesh <- trimesh +
  #ggtitle("(a)") +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme_light() +
  theme(legend.position = "none", plot.title = element_text(size = 5, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()#change legend key width
  )

trimesh_gr <- colour_long_edges(.data = distance, benchmark_value = benchmark,
                                triangular_object = tr1_object, distance_col = distance)

trimesh_gr <- trimesh_gr +
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

trimesh_removed <- remove_long_edges(.data = distance, benchmark_value = benchmark,
                                     triangular_object = tr1_object, distance_col = distance)

trimesh_removed <- trimesh_removed +
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



tour1 <- show_langevitour(df_all, df_bin, df_bin_centroids, benchmark_value = benchmark, distance = distance, distance_col = distance)

## To plot the distribution of distance
plot_dist <- function(distance_df){
  distance_df$group <- "1"
  dist_plot <- ggplot(distance_df, aes(x = group, y = distance)) +
    geom_quasirandom()+
    ylim(0, max(unlist(distance_df$distance))+ 0.5) + coord_flip()
  return(dist_plot)
}

distance_plot <- plot_dist(distance) +
  #ggtitle("(b)" ) +
  ylab(expression(d^{(2)})) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 12))

dist_p <- plot_dist(distance)  +
  geom_hline(yintercept = benchmark, colour = "red") +
  ylab(expression(d^{(2)})) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 12))

dist_p

trimesh + trimesh_gr +
  plot_layout(ncol = 2)

