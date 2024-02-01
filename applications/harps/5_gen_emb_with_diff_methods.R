library(dplyr)
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

set.seed(20240110)

use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

reticulate::source_python(paste0(here::here(), "/examples/function_scripts/Fit_PacMAP_code.py"))
reticulate::source_python(paste0(here::here(), "/examples/function_scripts/Fit_TriMAP_code.py"))

source("nldr_code.R", local = TRUE)
source("quollr_code.R", local = TRUE)

training_data <- read_rds("data/harps/harps_data.rds")
training_data <- training_data |>
  mutate(ID = 1:NROW(training_data))

clusters <- training_data$cluster
select_names <-  names(training_data)[-c(18, length(training_data))]

training_data <- training_data |> select(-cluster)
names(training_data) <- append(paste0("x", 1:17), "ID")

## tSNE (2)

tSNE_fit <- training_data |>
  select(-ID) |>
  Rtsne(perplexity = 91)

tSNE_data <- tSNE_fit$Y |>
  as.data.frame()

names(tSNE_data)[1:2] <- c("tSNE1", "tSNE2")

tSNE_data <- tSNE_data |>
  mutate(ID = training_data$ID)

## Run only once
write_rds(tSNE_data, "data/harps/harps_tsne_91.rds")

tSNE_data <- read_rds("data/harps/harps_tsne_91.rds")

plot_tSNE_2D(tSNE_data) + #ggtitle("(b)") +
  theme_linedraw() +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate(geom = 'text', label = 'b', x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 3)

num_bins_zeisel <- calculate_effective_x_bins(.data = tSNE_data, x = tSNE1,
                                              cell_area = 1) ##70

# num_bins_zeisel <- 30

shape_value_zeisel <- calculate_effective_shape_value(.data = tSNE_data,
                                                      x = tSNE1, y = tSNE2) ## 0.7831401

## To extract bin centroids
hexbin_data_object_zeisel <- extract_hexbin_centroids(nldr_df = tSNE_data,
                                                      num_bins = num_bins_zeisel,
                                                      shape_val = shape_value_zeisel, x = tSNE1, y = tSNE2)

df_bin_centroids_zeisel <- hexbin_data_object_zeisel$hexdf_data

##########


## Data set with all possible centroids in the hexagonal grid

full_centroid_df <- generate_full_grid_centroids(df_bin_centroids_zeisel)

## Generate all coordinates of hexagons
hex_grid <- full_hex_grid(full_centroid_df)

hex_full_count_df <- generate_full_grid_info(df_bin_centroids_zeisel)

#write_rds(hex_full_count_df, "data/zeisel/zeisel_tsne_hex_70.rds")

##########

ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_point(data = tSNE_data, aes(x = tSNE1, y = tSNE2), alpha = 0.5) +
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


# identify_rm_bins_zeisel <- find_low_density_hexagons(df_bin_centroids_zeisel, num_bins_zeisel)
#
# df_bin_centroids_zeisel <- df_bin_centroids_zeisel |>
#   filter(!(hexID %in% identify_rm_bins_zeisel))

tSNE_data_with_hb_id <- tSNE_data |>
  dplyr::mutate(hb_id = hexbin_data_object_zeisel$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_zeisel <- dplyr::bind_cols(training_data |> dplyr::select(-ID), tSNE_data_with_hb_id)

## Averaged on high-D
df_bin_zeisel <- avg_highD_data(.data = df_all_zeisel, column_start_text = "x")

## Triangulate bin centroids
tr1_object_zeisel <- triangulate_bin_centroids(df_bin_centroids_zeisel, x, y)
tr_from_to_df_zeisel <- generate_edge_info(triangular_object = tr1_object_zeisel)

## Compute 2D distances
distance_zeisel <- cal_2D_dist(.data = tr_from_to_df_zeisel)

plot_dist <- function(distance_df){
  distance_df$group <- "1"
  dist_plot <- ggplot(distance_df, aes(x = group, y = distance)) +
    geom_quasirandom()+
    ylim(0, max(unlist(distance_df$distance))+ 0.5) + coord_flip()
  return(dist_plot)
}

plot_dist(distance_zeisel) +
  #ggtitle("(b)" ) +
  ylab(expression(d^{(2)})) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 12))

## To find the benchmark value
benchmark_zeisel <- find_benchmark_value(.data = distance_zeisel, distance_col = distance)
#benchmark_zeisel <- 1.900017
##1.8

trimesh_zeisel_tsne <- ggplot(df_bin_centroids_zeisel, aes(x = x, y = y)) +
  geom_segment(data = tr_from_to_df_zeisel, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
  geom_point(size = 2, colour = "#33a02c") +
  coord_equal()

# ggplot(df_bin_centroids, aes(x = x, y = y)) +
# geom_point(size = 1, colour = "#33a02c") +
# geom_trimesh() +
# coord_equal()

trimesh_zeisel_tsne <- trimesh_zeisel_tsne +
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

trimesh_gr_zeisel_tsne <- colour_long_edges(.data = distance_zeisel, benchmark_value = benchmark_zeisel,
                                            triangular_object = tr1_object_zeisel, distance_col = distance)

trimesh_gr_zeisel_tsne <- trimesh_gr_zeisel_tsne +
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

trimesh_removed_zeisel_tsne <- remove_long_edges(.data = distance_zeisel, benchmark_value = benchmark_zeisel,
                                                 triangular_object = tr1_object_zeisel, distance_col = distance)

trimesh_removed_zeisel_tsne <- trimesh_removed_zeisel_tsne +
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

df_bin_train <- df_bin_zeisel
names(df_bin_train)[-1] <- paste0("avg_", names(df_bin_train)[-1])

error_df <- df_all_zeisel |>
  dplyr::left_join(df_bin_train, by = c("hb_id" = "hb_id")) ## Map high-D averaged/weighted mean coordinates

# prediction_df_join <- prediction_df_join |>
#   dplyr::left_join(data, by = c("ID" = "ID")) ## Map high-D data

for (i in 1:(NCOL(df_bin_train) - 1)) {

  error_df[ , paste0("abs_residual_", "x", i)] <- abs(error_df[ , paste0("x", i)] - error_df[ , paste0("avg_", "x", i)])

}

error_df <- error_df |>
  dplyr::mutate(total = rowSums(dplyr::pick(tidyselect::starts_with(paste0("abs_residual_", "x")))))

library(ggbeeswarm)
error_df$group <- "1"
ggplot(error_df, aes(x = group, y = total)) +
  geom_quasirandom()+
  ylim(0, max(unlist(error_df$total))+ 0.5) + coord_flip()

### The minimum error is 0 and the maximum is 1.38
### There is lot of points with error 0,

error_df <- error_df |>
  mutate(type = if_else(total <= 0, "no error",
                        if_else(total <= 0.3, "error 0-0.3",
                                if_else(total <= 0.5, "error 0.3-0.5",
                                        if_else(total <= 1, "error 0.5-1",
                                                if_else(total <= 1.5, "error 1-1.5",
                                                        "error greter than 1.5"))))))


####
### Define type column
df <- df_all_zeisel |>
  dplyr::select(tidyselect::starts_with("x")) #|>
#dplyr::mutate(type = "data")  ## original dataset

residual_df <- error_df |> select(type)

df <- dplyr::bind_cols(df, residual_df)

df_b <- df_bin_zeisel |>
  dplyr::filter(hb_id %in% df_bin_centroids_zeisel$hexID) |>
  dplyr::select(-hb_id) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

df_exe <- dplyr::bind_rows(df_b, df)
names(df_exe) <- append(select_names, "type")

distance_df_small_edges <- distance_zeisel %>%
  dplyr::filter(distance < benchmark_zeisel)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to, group = factor(df_exe$type, levels = c("no error", "error 0-0.3", "error 0.3-0.5", "error 0.5-1", "error 1-1.5", "error greter than 1.5", "model")), pointSize = 3,
                         levelColors = c("#b15928", "#1f78b4",
                                         "#fb9a99", "#cab2d6", "#ff7f00",
                                         "#e31a1c", "#33a02c"))



#tour_umap_harps <- show_langevitour(df_all_umap_harps, df_bin_umap_harps, df_bin_centroids_umap_harps, benchmark_value = benchmark_umap_harps, distance = distance_umap_harps, distance_col = distance)


## tSNE (1)

tSNE_fit <- training_data |>
  select(-ID) |>
  Rtsne(perplexity = 23)

tSNE_data <- tSNE_fit$Y |>
  as.data.frame()

names(tSNE_data)[1:2] <- c("tSNE1", "tSNE2")

tSNE_data <- tSNE_data |>
  mutate(ID = training_data$ID)

## Run only once
write_rds(tSNE_data, "data/harps/harps_tsne_23.rds")

tSNE_data <- read_rds("data/harps/harps_tsne_23.rds")

plot_tSNE_2D(tSNE_data) + #ggtitle("(b)") +
  theme_linedraw() +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate(geom = 'text', label = 'b', x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 3)

num_bins_zeisel <- calculate_effective_x_bins(.data = tSNE_data, x = tSNE1,
                                              cell_area = 1) ##70

num_bins_zeisel <- 30

shape_value_zeisel <- calculate_effective_shape_value(.data = tSNE_data,
                                                      x = tSNE1, y = tSNE2) ## 0.7831401

## To extract bin centroids
hexbin_data_object_zeisel <- extract_hexbin_centroids(nldr_df = tSNE_data,
                                                      num_bins = num_bins_zeisel,
                                                      shape_val = shape_value_zeisel, x = tSNE1, y = tSNE2)

df_bin_centroids_zeisel <- hexbin_data_object_zeisel$hexdf_data

##########


## Data set with all possible centroids in the hexagonal grid

full_centroid_df <- generate_full_grid_centroids(df_bin_centroids_zeisel)

## Generate all coordinates of hexagons
hex_grid <- full_hex_grid(full_centroid_df)

hex_full_count_df <- generate_full_grid_info(df_bin_centroids_zeisel)

#write_rds(hex_full_count_df, "data/zeisel/zeisel_tsne_hex_70.rds")

##########

ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_point(data = tSNE_data, aes(x = tSNE1, y = tSNE2), alpha = 0.5) +
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


# identify_rm_bins_zeisel <- find_low_density_hexagons(df_bin_centroids_zeisel, num_bins_zeisel)
#
# df_bin_centroids_zeisel <- df_bin_centroids_zeisel |>
#   filter(!(hexID %in% identify_rm_bins_zeisel))

tSNE_data_with_hb_id <- tSNE_data |>
  dplyr::mutate(hb_id = hexbin_data_object_zeisel$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_zeisel <- dplyr::bind_cols(training_data |> dplyr::select(-ID), tSNE_data_with_hb_id)

## Averaged on high-D
df_bin_zeisel <- avg_highD_data(.data = df_all_zeisel, column_start_text = "x")

## Triangulate bin centroids
tr1_object_zeisel <- triangulate_bin_centroids(df_bin_centroids_zeisel, x, y)
tr_from_to_df_zeisel <- generate_edge_info(triangular_object = tr1_object_zeisel)

## Compute 2D distances
distance_zeisel <- cal_2D_dist(.data = tr_from_to_df_zeisel)

plot_dist <- function(distance_df){
  distance_df$group <- "1"
  dist_plot <- ggplot(distance_df, aes(x = group, y = distance)) +
    geom_quasirandom()+
    ylim(0, max(unlist(distance_df$distance))+ 0.5) + coord_flip()
  return(dist_plot)
}

plot_dist(distance_zeisel) +
  #ggtitle("(b)" ) +
  ylab(expression(d^{(2)})) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 12))

## To find the benchmark value
benchmark_zeisel <- find_benchmark_value(.data = distance_zeisel, distance_col = distance)
#benchmark_zeisel <- 1.900017
##1.8

trimesh_zeisel_tsne <- ggplot(df_bin_centroids_zeisel, aes(x = x, y = y)) +
  geom_segment(data = tr_from_to_df_zeisel, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
  geom_point(size = 2, colour = "#33a02c") +
  coord_equal()

# ggplot(df_bin_centroids, aes(x = x, y = y)) +
# geom_point(size = 1, colour = "#33a02c") +
# geom_trimesh() +
# coord_equal()

trimesh_zeisel_tsne <- trimesh_zeisel_tsne +
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

trimesh_gr_zeisel_tsne <- colour_long_edges(.data = distance_zeisel, benchmark_value = benchmark_zeisel,
                                            triangular_object = tr1_object_zeisel, distance_col = distance)

trimesh_gr_zeisel_tsne <- trimesh_gr_zeisel_tsne +
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

trimesh_removed_zeisel_tsne <- remove_long_edges(.data = distance_zeisel, benchmark_value = benchmark_zeisel,
                                                 triangular_object = tr1_object_zeisel, distance_col = distance)

trimesh_removed_zeisel_tsne <- trimesh_removed_zeisel_tsne +
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

df_bin_train <- df_bin_zeisel
names(df_bin_train)[-1] <- paste0("avg_", names(df_bin_train)[-1])

error_df <- df_all_zeisel |>
  dplyr::left_join(df_bin_train, by = c("hb_id" = "hb_id")) ## Map high-D averaged/weighted mean coordinates

# prediction_df_join <- prediction_df_join |>
#   dplyr::left_join(data, by = c("ID" = "ID")) ## Map high-D data

for (i in 1:(NCOL(df_bin_train) - 1)) {

  error_df[ , paste0("abs_residual_", "x", i)] <- abs(error_df[ , paste0("x", i)] - error_df[ , paste0("avg_", "x", i)])

}

error_df <- error_df |>
  dplyr::mutate(total = rowSums(dplyr::pick(tidyselect::starts_with(paste0("abs_residual_", "x")))))

library(ggbeeswarm)
error_df$group <- "1"
ggplot(error_df, aes(x = group, y = total)) +
  geom_quasirandom()+
  ylim(0, max(unlist(error_df$total))+ 0.5) + coord_flip()

### The minimum error is 0 and the maximum is 1.38
### There is lot of points with error 0,

error_df <- error_df |>
  mutate(type = if_else(total <= 0, "no error",
                        if_else(total <= 0.125, "error 0-0.125",
                                if_else(total <= 0.25, "error 0.125-0.25",
                                        if_else(total <= 0.5, "error 0.25-0.5",
                                                if_else(total <= 1, "error 0.5-1",
                                                        "error greter than 1"))))))


####
### Define type column
df <- df_all_zeisel |>
  dplyr::select(tidyselect::starts_with("x")) #|>
#dplyr::mutate(type = "data")  ## original dataset

residual_df <- error_df |> select(type)

df <- dplyr::bind_cols(df, residual_df)

df_b <- df_bin_zeisel |>
  dplyr::filter(hb_id %in% df_bin_centroids_zeisel$hexID) |>
  dplyr::select(-hb_id) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

df_exe <- dplyr::bind_rows(df_b, df)
names(df_exe) <- append(select_names, "type")

distance_df_small_edges <- distance_zeisel %>%
  dplyr::filter(distance < benchmark_zeisel)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to, group = factor(df_exe$type, levels = c("no error", "error 0-0.125", "error 0.125-0.25", "error 0.25-0.5", "error 0.5-1", "error greter than 1", "model")), pointSize = 3,
                         levelColors = c("#b15928", "#1f78b4",
                                         "#fb9a99", "#cab2d6", "#ff7f00",
                                         "#e31a1c", "#33a02c"))



#tour_umap_harps <- show_langevitour(df_all_umap_harps, df_bin_umap_harps, df_bin_centroids_umap_harps, benchmark_value = benchmark_umap_harps, distance = distance_umap_harps, distance_col = distance)



## UMAP

UMAP_fit <- umap(training_data |> dplyr::select(-ID), n_neighbors = 15, n_components =  2, min_dist = 0.001)

UMAP_data <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data)[1:(ncol(UMAP_data))] <- paste0(rep("UMAP",(ncol(UMAP_data))), 1:(ncol(UMAP_data)))

UMAP_data <- UMAP_data |>
  mutate(ID = training_data$ID)

plot_UMAP_2D(UMAP_data) + #ggtitle("(b)") +
  theme_linedraw() +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate(geom = 'text', label = 'b', x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 3)

num_bins_pbmc <- calculate_effective_x_bins(.data = UMAP_data, x = UMAP1,
                                            cell_area = 1) ##23

shape_value_pbmc <- calculate_effective_shape_value(.data = UMAP_data,
                                                    x = UMAP1, y = UMAP2) ## 0.8772751

## To extract bin centroids
hexbin_data_object_pbmc <- extract_hexbin_centroids(nldr_df = UMAP_data,
                                                    num_bins = num_bins_pbmc,
                                                    shape_val = shape_value_pbmc, x = UMAP1, y = UMAP2)

df_bin_centroids_pbmc <- hexbin_data_object_pbmc$hexdf_data

##########

## Data set with all possible centroids in the hexagonal grid

full_centroid_df <- generate_full_grid_centroids(df_bin_centroids_pbmc)

## Generate all coordinates of hexagons
hex_grid <- full_hex_grid(full_centroid_df)

hex_full_count_df <- generate_full_grid_info(df_bin_centroids_pbmc)

#write_rds(hex_full_count_df, "data/pbmc/pbmc_3k_festem/pbmc_umap_hex_23.rds")

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


UMAP_pbmc_with_hb_id <- UMAP_data |>
  dplyr::mutate(hb_id = hexbin_data_object_pbmc$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_pbmc <- dplyr::bind_cols(training_data|> dplyr::select(-ID), UMAP_pbmc_with_hb_id)

## Averaged on high-D
df_bin_pbmc <- avg_highD_data(.data = df_all_pbmc, column_start_text = "x")

## Triangulate bin centroids
tr1_object_pbmc <- triangulate_bin_centroids(df_bin_centroids_pbmc, x, y)
tr_from_to_df_pbmc <- generate_edge_info(triangular_object = tr1_object_pbmc)

## Compute 2D distances
distance_pbmc <- cal_2D_dist(.data = tr_from_to_df_pbmc)

## To find the benchmark value
benchmark_pbmc <- find_benchmark_value(.data = distance_pbmc, distance_col = distance)
benchmark_pbmc <- 0.9774252


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

df_bin_train <- df_bin_pbmc
names(df_bin_train)[-1] <- paste0("avg_", names(df_bin_train)[-1])

error_df <- df_all_pbmc |>
  dplyr::left_join(df_bin_train, by = c("hb_id" = "hb_id")) ## Map high-D averaged/weighted mean coordinates

# prediction_df_join <- prediction_df_join |>
#   dplyr::left_join(data, by = c("ID" = "ID")) ## Map high-D data

for (i in 1:(NCOL(df_bin_train) - 1)) {

  error_df[ , paste0("abs_residual_", "x", i)] <- abs(error_df[ , paste0("x", i)] - error_df[ , paste0("avg_", "x", i)])

}

error_df <- error_df |>
  dplyr::mutate(total = rowSums(dplyr::pick(tidyselect::starts_with(paste0("abs_residual_", "x")))))

library(ggbeeswarm)
error_df$group <- "1"
ggplot(error_df, aes(x = group, y = total)) +
  geom_quasirandom()+
  ylim(0, max(unlist(error_df$total))+ 0.5) + coord_flip()

### The minimum error is 0 and the maximum is 42.17439
### There is lot of points with error 0,

error_df <- error_df |>
  mutate(type = if_else(total <= 0, "no error",
                        if_else(total <= 1, "error 0-1",
                                if_else(total <= 2, "error 1-2",
                                        if_else(total <= 3, "error 2-3",
                                                if_else(total <= 4, "error 3-4",
                                                        "error greter than 4"))))))

### Define type column
df <- df_all_pbmc |>
  dplyr::select(tidyselect::starts_with("x")) #|>
#dplyr::rename("type" = "cell_label") ## original dataset

residual_df <- error_df |> select(type)

df <- dplyr::bind_cols(df, residual_df)


#df$type <- as.factor(df$type)

#levels(df$type) <- c("Memory \nCD4 T", "Naive CD4 T", "CD14+ Mono",  "B", "CD8 T", "FCGR3A+ \n Mono", "NK", "M-MDSC\n-like", "CD27-CD+ \n Memory T", "DC")

df_b <- df_bin_pbmc |>
  dplyr::filter(hb_id %in% df_bin_centroids_pbmc$hexID) |>
  dplyr::select(-hb_id) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

df_exe <- dplyr::bind_rows(df_b, df)

distance_df_small_edges <- distance_pbmc %>%
  dplyr::filter(distance < benchmark_pbmc)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to, group = factor(df_exe$type , levels = c("no error", "error 0-1", "error 1-2", "error 2-3", "error 3-4", "error greter than 4", "model")
                         ), pointSize = 3,
                         levelColors = c("#b15928", "#1f78b4",
                                         "#fb9a99", "#cab2d6", "#ff7f00",
                                         "#e31a1c", "#33a02c"))


tour1_pbmc_umap <- show_langevitour(df_all_pbmc, df_bin_pbmc,
                                    df_bin_centroids_pbmc, benchmark_value = benchmark_pbmc,
                                    distance = distance_pbmc, distance_col = distance, col_start = "x")



## Run only once
#write_rds(UMAP_data, file = "data/s_curve/s_curve_umap.rds")

# predict_UMAP_df <- predict(UMAP_fit, test_data |> dplyr::select(-ID)) |>
#   as.data.frame()
#
# names(predict_UMAP_df)[1:(ncol(predict_UMAP_df))] <- paste0(rep("UMAP",(ncol(predict_UMAP_df))), 1:(ncol(predict_UMAP_df)))
#
# predict_UMAP_df <- predict_UMAP_df |>
#   mutate(ID = test_data$ID)

## Run only once
#write_rds(UMAP_data, file = "data/s_curve/s_curve_umap_predict.rds")

## PHATE

PHATE_data <- Fit_PHATE(training_data |> dplyr::select(-ID), knn = 5, with_seed = 20240110)
PHATE_data <- PHATE_data |>
  select(PHATE1, PHATE2)
PHATE_data <- PHATE_data |>
  mutate(ID = training_data$ID)

#write_rds(PHATE_data, file = "data/s_curve/s_curve_phate.rds")

plot_PHATE_2D(PHATE_data) + #ggtitle("(c)") +
  theme_linedraw() +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate(geom = 'text', label = 'c', x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 3)


## TriMAP

tem_dir <- tempdir()

Fit_TriMAP_data(training_data |> dplyr::select(-ID), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_TriMAP_values.csv")

Fit_TriMAP(as.integer(2), as.integer(5), as.integer(4), as.integer(3), path, path2)

TriMAP_data <- read_csv(path2)
TriMAP_data <- TriMAP_data |>
  mutate(ID = training_data$ID)
#write_rds(TriMAP_data, file = "data/s_curve/s_curve_trimap.rds")

plot_TriMAP_2D(TriMAP_data) + #ggtitle("(d)") +
  theme_linedraw() +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate(geom = 'text', label = 'd', x = Inf, y = Inf, hjust = 1.5, vjust = 1.5, size = 3)

num_bins_pbmc <- calculate_effective_x_bins(.data = TriMAP_data, x = TriMAP1,
                                            cell_area = 1) ##23

shape_value_pbmc <- calculate_effective_shape_value(.data = TriMAP_data,
                                                    x = TriMAP1, y = TriMAP2) ## 0.8772751

## To extract bin centroids
hexbin_data_object_pbmc <- extract_hexbin_centroids(nldr_df = TriMAP_data,
                                                    num_bins = num_bins_pbmc,
                                                    shape_val = shape_value_pbmc, x = TriMAP1, y = TriMAP2)

df_bin_centroids_pbmc <- hexbin_data_object_pbmc$hexdf_data

##########

## Data set with all possible centroids in the hexagonal grid

full_centroid_df <- generate_full_grid_centroids(df_bin_centroids_pbmc)

## Generate all coordinates of hexagons
hex_grid <- full_hex_grid(full_centroid_df)

hex_full_count_df <- generate_full_grid_info(df_bin_centroids_pbmc)

#write_rds(hex_full_count_df, "data/pbmc/pbmc_3k_festem/pbmc_umap_hex_23.rds")

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


UMAP_pbmc_with_hb_id <- TriMAP_data |>
  dplyr::mutate(hb_id = hexbin_data_object_pbmc$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_pbmc <- dplyr::bind_cols(training_data|> dplyr::select(-ID), UMAP_pbmc_with_hb_id)

## Averaged on high-D
df_bin_pbmc <- avg_highD_data(.data = df_all_pbmc, column_start_text = "x")

## Triangulate bin centroids
tr1_object_pbmc <- triangulate_bin_centroids(df_bin_centroids_pbmc, x, y)
tr_from_to_df_pbmc <- generate_edge_info(triangular_object = tr1_object_pbmc)

## Compute 2D distances
distance_pbmc <- cal_2D_dist(.data = tr_from_to_df_pbmc)

## To find the benchmark value
benchmark_pbmc <- find_benchmark_value(.data = distance_pbmc, distance_col = distance)
benchmark_pbmc <- 2.131584


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

df_bin_train <- df_bin_pbmc
names(df_bin_train)[-1] <- paste0("avg_", names(df_bin_train)[-1])

error_df <- df_all_pbmc |>
  dplyr::left_join(df_bin_train, by = c("hb_id" = "hb_id")) ## Map high-D averaged/weighted mean coordinates

# prediction_df_join <- prediction_df_join |>
#   dplyr::left_join(data, by = c("ID" = "ID")) ## Map high-D data

for (i in 1:(NCOL(df_bin_train) - 1)) {

  error_df[ , paste0("abs_residual_", "x", i)] <- abs(error_df[ , paste0("x", i)] - error_df[ , paste0("avg_", "x", i)])

}

error_df <- error_df |>
  dplyr::mutate(total = rowSums(dplyr::pick(tidyselect::starts_with(paste0("abs_residual_", "x")))))

library(ggbeeswarm)
error_df$group <- "1"
ggplot(error_df, aes(x = group, y = total)) +
  geom_quasirandom()+
  ylim(0, max(unlist(error_df$total))+ 0.5) + coord_flip()

### The minimum error is 0 and the maximum is 42.17439
### There is lot of points with error 0,

error_df <- error_df |>
  mutate(type = if_else(total <= 0, "no error",
                        if_else(total <= 1, "error 0-1",
                                if_else(total <= 2, "error 1-2",
                                        if_else(total <= 3, "error 2-3",
                                                if_else(total <= 4, "error 3-4",
                                                        "error greter than 4"))))))

### Define type column
df <- df_all_pbmc |>
  dplyr::select(tidyselect::starts_with("x")) #|>
#dplyr::rename("type" = "cell_label") ## original dataset

residual_df <- error_df |> select(type)

df <- dplyr::bind_cols(df, residual_df)


#df$type <- as.factor(df$type)

#levels(df$type) <- c("Memory \nCD4 T", "Naive CD4 T", "CD14+ Mono",  "B", "CD8 T", "FCGR3A+ \n Mono", "NK", "M-MDSC\n-like", "CD27-CD+ \n Memory T", "DC")

df_b <- df_bin_pbmc |>
  dplyr::filter(hb_id %in% df_bin_centroids_pbmc$hexID) |>
  dplyr::select(-hb_id) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

df_exe <- dplyr::bind_rows(df_b, df)

distance_df_small_edges <- distance_pbmc %>%
  dplyr::filter(distance < benchmark_pbmc)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to, group = factor(df_exe$type , levels = c("no error", "error 0-1", "error 1-2", "error 2-3", "error 3-4", "error greter than 4", "model")
                         ), pointSize = 3,
                         levelColors = c("#b15928", "#1f78b4",
                                         "#fb9a99", "#cab2d6", "#ff7f00",
                                         "#e31a1c", "#33a02c"))


tour1_pbmc_umap <- show_langevitour(df_all_pbmc, df_bin_pbmc,
                                    df_bin_centroids_pbmc, benchmark_value = benchmark_pbmc,
                                    distance = distance_pbmc, distance_col = distance, col_start = "x")





## PaCMAP


tem_dir <- tempdir()

Fit_PacMAP_data(training_data |> dplyr::select(-ID), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_PaCMAP_values.csv")

Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.9, as.integer(2), path, path2)

PaCMAP_data <- read_csv(path2)
PaCMAP_data <- PaCMAP_data |>
  mutate(ID = training_data$ID)

#write_rds(PacMAP_data, file = "data/s_curve/s_curve_pacmap.rds")
PaCMAP_data |>
  ggplot(aes(x = PaCMAP1,
             y = PaCMAP2, colour = clusters))+
  geom_point(alpha=0.5) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + #ggtitle("(a)") +
  theme_linedraw()  +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm')) +
  scale_color_manual(values=c("#b15928", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#6a3d9a", "#ff7f00", "#cab2d6", "#fdbf6f", "#ffff99", "#a6cee3", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#000000", "#bdbdbd")) +
  guides(color = guide_legend(title = "cluster"),
         color = guide_legend(nrow = 3)) +
  annotate(geom = 'text', label = 'e', x = Inf, y = Inf, hjust = 1.5, vjust = 1.5, size = 3)

num_bins_pbmc <- calculate_effective_x_bins(.data = PaCMAP_data, x = PaCMAP1,
                                            cell_area = 1) ##23

shape_value_pbmc <- calculate_effective_shape_value(.data = PaCMAP_data,
                                                    x = PaCMAP1, y = PaCMAP2) ## 0.8772751

## To extract bin centroids
hexbin_data_object_pbmc <- extract_hexbin_centroids(nldr_df = PaCMAP_data,
                                                    num_bins = num_bins_pbmc,
                                                    shape_val = shape_value_pbmc, x = PaCMAP1, y = PaCMAP2)

df_bin_centroids_pbmc <- hexbin_data_object_pbmc$hexdf_data

##########

## Data set with all possible centroids in the hexagonal grid

full_centroid_df <- generate_full_grid_centroids(df_bin_centroids_pbmc)

## Generate all coordinates of hexagons
hex_grid <- full_hex_grid(full_centroid_df)

hex_full_count_df <- generate_full_grid_info(df_bin_centroids_pbmc)

#write_rds(hex_full_count_df, "data/pbmc/pbmc_3k_festem/pbmc_umap_hex_23.rds")

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


UMAP_pbmc_with_hb_id <- PaCMAP_data |>
  dplyr::mutate(hb_id = hexbin_data_object_pbmc$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_pbmc <- dplyr::bind_cols(training_data|> dplyr::select(-ID), UMAP_pbmc_with_hb_id)

## Averaged on high-D
df_bin_pbmc <- avg_highD_data(.data = df_all_pbmc, column_start_text = "x")

## Triangulate bin centroids
tr1_object_pbmc <- triangulate_bin_centroids(df_bin_centroids_pbmc, x, y)
tr_from_to_df_pbmc <- generate_edge_info(triangular_object = tr1_object_pbmc)

## Compute 2D distances
distance_pbmc <- cal_2D_dist(.data = tr_from_to_df_pbmc)

## To find the benchmark value
benchmark_pbmc <- find_benchmark_value(.data = distance_pbmc, distance_col = distance)
benchmark_pbmc <- 1.020554


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

df_bin_train <- df_bin_pbmc
names(df_bin_train)[-1] <- paste0("avg_", names(df_bin_train)[-1])

error_df <- df_all_pbmc |>
  dplyr::left_join(df_bin_train, by = c("hb_id" = "hb_id")) ## Map high-D averaged/weighted mean coordinates

# prediction_df_join <- prediction_df_join |>
#   dplyr::left_join(data, by = c("ID" = "ID")) ## Map high-D data

for (i in 1:(NCOL(df_bin_train) - 1)) {

  error_df[ , paste0("abs_residual_", "x", i)] <- abs(error_df[ , paste0("x", i)] - error_df[ , paste0("avg_", "x", i)])

}

error_df <- error_df |>
  dplyr::mutate(total = rowSums(dplyr::pick(tidyselect::starts_with(paste0("abs_residual_", "x")))))

library(ggbeeswarm)
error_df$group <- "1"
ggplot(error_df, aes(x = group, y = total)) +
  geom_quasirandom()+
  ylim(0, max(unlist(error_df$total))+ 0.5) + coord_flip()

### The minimum error is 0 and the maximum is 42.17439
### There is lot of points with error 0,

error_df <- error_df |>
  mutate(type = if_else(total <= 0, "no error",
                        if_else(total <= 1, "error 0-1",
                                if_else(total <= 2, "error 1-2",
                                        if_else(total <= 3, "error 2-3",
                                                if_else(total <= 4, "error 3-4",
                                                        "error greter than 4"))))))

### Define type column
df <- df_all_pbmc |>
  dplyr::select(tidyselect::starts_with("x")) #|>
#dplyr::rename("type" = "cell_label") ## original dataset

residual_df <- error_df |> select(type)

df <- dplyr::bind_cols(df, residual_df)


#df$type <- as.factor(df$type)

#levels(df$type) <- c("Memory \nCD4 T", "Naive CD4 T", "CD14+ Mono",  "B", "CD8 T", "FCGR3A+ \n Mono", "NK", "M-MDSC\n-like", "CD27-CD+ \n Memory T", "DC")

df_b <- df_bin_pbmc |>
  dplyr::filter(hb_id %in% df_bin_centroids_pbmc$hexID) |>
  dplyr::select(-hb_id) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

df_exe <- dplyr::bind_rows(df_b, df)

distance_df_small_edges <- distance_pbmc %>%
  dplyr::filter(distance < benchmark_pbmc)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to, group = factor(df_exe$type , levels = c("no error", "error 0-1", "error 1-2", "error 2-3", "error 3-4", "error greter than 4", "model")
                         ), pointSize = 3,
                         levelColors = c("#b15928", "#1f78b4",
                                         "#fb9a99", "#cab2d6", "#ff7f00",
                                         "#e31a1c", "#33a02c"))


tour1_pbmc_umap <- show_langevitour(df_all_pbmc, df_bin_pbmc,
                                    df_bin_centroids_pbmc, benchmark_value = benchmark_pbmc,
                                    distance = distance_pbmc, distance_col = distance, col_start = "x")



