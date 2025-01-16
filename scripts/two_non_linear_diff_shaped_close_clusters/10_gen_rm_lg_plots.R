## This script is to generate tuning the long edge with UMAP layout

## Import necessary libraries
library(quollr)
library(dplyr)
library(readr)
library(langevitour)
library(conflicted)
library(ggplot2)
library(colorspace)

conflicts_prefer(dplyr::filter)


training_data_two_curvy <- read_rds("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_data.rds")
umap_two_curvy <- read_rds(file = "data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_umap_n-neigbors_15_min-dist_0.1.rds")

data_two_curvy <- training_data_two_curvy |>
  mutate(type = "data")

two_curvy_scaled_obj <- gen_scaled_data(
  data = umap_two_curvy)

umap_two_curvy_scaled <- two_curvy_scaled_obj$scaled_nldr |>
  mutate(ID = row_number())
lim1 <- two_curvy_scaled_obj$lim1
lim2 <- two_curvy_scaled_obj$lim2
r2 <- diff(lim2)/diff(lim1)

sc_ltr_pos <- c(0.08, 0.9)

##Full hexagon grid with UMAP data

## Compute hexbin parameters
num_bins_x_two_curvy <- 42

## hexagon binning to have regular hexagons
hb_obj_two_curvy <- hex_binning(
  data = umap_two_curvy_scaled,
  bin1 = num_bins_x_two_curvy,
  r2 = r2,
  q = 0.1)

## Data set with all centroids
all_centroids_df <- hb_obj_two_curvy$centroids

a1_1 <- hb_obj_two_curvy$a1

## Generate all coordinates of hexagons
hex_grid <- hb_obj_two_curvy$hex_poly

## To obtain the standardise counts within hexbins
counts_df <- hb_obj_two_curvy$std_cts
df_bin_centroids_two_curvy <- extract_hexbin_centroids(
  centroids_df = all_centroids_df,
  counts_df = counts_df) |>
  filter(drop_empty == FALSE)

umap_data_two_curvy_with_hb_id <- hb_obj_two_curvy$data_hb_id
df_all_two_curvy <- dplyr::bind_cols(training_data_two_curvy, umap_data_two_curvy_with_hb_id)
df_bin_two_curvy <- avg_highd_data(data = df_all_two_curvy, col_start = "x")

hex_grid_with_counts <-
  left_join(hex_grid,
            counts_df,
            by = c("hex_poly_id" = "hb_id"))

##Non-empty bins with bin centroids

hex_grid_nonempty <- hex_grid |>
  filter(hex_poly_id %in% df_bin_centroids_two_curvy$hexID)

##2D model

## Triangulate bin centroids
tr1_object_two_curvy <- tri_bin_centroids(
  df_bin_centroids_two_curvy, x = "c_x", y = "c_y")
tr_from_to_df_two_curvy <- gen_edges(
  tri_object = tr1_object_two_curvy)

trimesh_two_curvy_umap <- ggplot() +
  geom_segment(data = tr_from_to_df_two_curvy,
               aes(
                 x = x_from,
                 y = y_from,
                 xend = x_to,
                 yend = y_to),
               colour = "#000000") +
  geom_point(data = umap_two_curvy_scaled,
             aes(
               x = UMAP1,
               y = UMAP2
             ),
             color = "#636363",
             alpha = 0.5,
             size = 0.5
  ) +
  coord_equal() +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  interior_annotation("a1", sc_ltr_pos)

## Compute 2D distances
distance_two_curvy <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_two_curvy,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_two_curvy <- find_lg_benchmark(
  distance_edges = distance_two_curvy,
  distance_col = "distance")

distance_df_small_edges_two_curvy <- distance_two_curvy |>
  filter(distance < benchmark_two_curvy)

bin_width <- hb_obj_two_curvy$a1

####With diff benchmark values

## Benchmark 2

benchmark_two_curvy <- 6.5 * bin_width

distance_df_small_edges_two_curvy3 <- distance_two_curvy |>
  filter(distance < benchmark_two_curvy)

tr_from_to_df_two_curvy4 <- inner_join(
  tr_from_to_df_two_curvy, distance_df_small_edges_two_curvy3,
  by = c("from", "to"))

trimesh_two_curvy_removed2_umap <- ggplot() +
  geom_segment(data = tr_from_to_df_two_curvy4,
               aes(
                 x = x_from,
                 y = y_from,
                 xend = x_to,
                 yend = y_to),
               colour = "#000000")+
  geom_point(data = umap_two_curvy_scaled,
             aes(
               x = UMAP1,
               y = UMAP2
             ),
             color = "#636363",
             alpha = 0.5,
             size = 0.5
  ) +
  coord_equal() +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  interior_annotation("b1", sc_ltr_pos)

## Benchmark 1

benchmark_two_curvy <- 2 * bin_width

distance_df_small_edges_two_curvy2 <- distance_two_curvy |>
  filter(distance < benchmark_two_curvy) ## 0.231

tr_from_to_df_two_curvy3 <- inner_join(
  tr_from_to_df_two_curvy, distance_df_small_edges_two_curvy2,
  by = c("from", "to"))

trimesh_two_curvy_removed1_umap <- ggplot() +
  geom_segment(data = tr_from_to_df_two_curvy3,
               aes(
                 x = x_from,
                 y = y_from,
                 xend = x_to,
                 yend = y_to),
               colour = "#000000")+
  geom_point(data = umap_two_curvy_scaled,
             aes(
               x = UMAP1,
               y = UMAP2
             ),
             color = "#636363",
             alpha = 0.5,
             size = 0.5
  ) +
  coord_equal() +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  interior_annotation("c1", sc_ltr_pos)

trimesh_two_curvy_removed1_with_data <- ggplot() +
  geom_segment(data = tr_from_to_df_two_curvy3,
               aes(
                 x = x_from,
                 y = y_from,
                 xend = x_to,
                 yend = y_to),
               colour = "#000000") +
  geom_point(data = umap_two_curvy_scaled,
             aes(
               x = UMAP1,
               y = UMAP2
             ),
             color = "#636363",
             alpha = 0.5,
             size = 0.5
  ) +
  coord_equal() +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  interior_annotation("a2", sc_ltr_pos)

## Computed benchmark value
tr_from_to_df_two_curvy <- inner_join(
  tr_from_to_df_two_curvy, distance_df_small_edges_two_curvy,
  by = c("from", "to"))

df_bin_two_curvy_temp <- df_bin_two_curvy

### Compute highD distance
dist_vec <- proxy::dist(x = df_bin_two_curvy_temp[, -1], method = "Euclidean") |> as.vector()

from_vec <- c()
to_vec <- c()
num_obs <- 1:(NROW(df_bin_two_curvy_temp) - 1)

for (obs in num_obs) {

  from_val <- rep(obs, (NROW(df_bin_two_curvy_temp) - obs))
  if ((obs + 1) <= NROW(df_bin_two_curvy_temp)) {
    to_val <- (obs + 1):NROW(df_bin_two_curvy_temp)
  }
  from_vec <- append(from_vec, from_val)
  to_vec <- append(to_vec, to_val)

}

dist_highd <- tibble::tibble(from = from_vec, to = to_vec, dist_highd = dist_vec)


## To plot the distribution of distances

# Define bin width and starting point
start_point <- bin_width/2

# Create bins and calculate the mean value for each bin range
bins <- seq(start_point, max(distance_two_curvy$distance) + bin_width,
            by = bin_width)
bin_labels <- bins[-length(bins)] + bin_width / 2  # mean of each bin range

distance_two_curvy_dist <- distance_two_curvy |>
  mutate(bin = cut(distance, breaks = bins, include.lowest = TRUE, labels = bin_labels))

distance_two_curvy_dist <- left_join(distance_two_curvy_dist, dist_highd, by = c("from", "to"))

text_df <- tibble::tibble(
  x = c(0.035, 0.15),
  y = c(10, 10),
  text = c("c", "b")
)

distance_points_umap <- ggplot(
  distance_two_curvy_dist,
  aes(x = distance,
      y = dist_highd)) +
  geom_point(alpha = 0.5) +
  geom_text(
    data=text_df,
    aes(x=x, y=y,
        label=text),
    colour="#bdbdbd",
    size = 10) +
  geom_vline(xintercept = 6.5 * bin_width,
             linetype="solid",
             color = "#bdbdbd",
             linewidth=1) +
  geom_vline(xintercept = 2 * bin_width,
             linetype="solid",
             color = "#bdbdbd",
             linewidth=1) +
  ylab(expression(d^{(7)})) +
  xlab(expression(d^{(2)})) +
  theme_minimal() +
  theme(aspect.ratio = 0.75,
        panel.border = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 12, hjust = 0.5, vjust = -0.5),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_line())


df_bin_two_curvy <- df_bin_two_curvy |>
  select(-hb_id) |>
  mutate(type = "model")

# Apply the scaling
df_model_data_two_curvy <- bind_rows(data_two_curvy, df_bin_two_curvy)

scaled_two_curvy <- scale_data_manual(df_model_data_two_curvy, "type") |>
  as_tibble()

scaled_two_curvy_data <- scaled_two_curvy |>
  filter(type == "data") |>
  select(-type)

scaled_two_curvy_data_model <- scaled_two_curvy |>
  filter(type == "model") |>
  select(-type)


## First projection
projection <- cbind(
  c(-0.02291,-0.00834,0.06450,-0.09921,0.00761,-0.05518,-0.01158),
  c(-0.11183,-0.02761,-0.04394,-0.01066,0.04229,0.02054,0.01768))

projection_scaled <- projection * 5

projected <- as.matrix(scaled_two_curvy_data) %*% projection_scaled

projected_df <- projected |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::rename(c("proj1" = "...1",
                  "proj2" = "...2")) |>
  dplyr::mutate(ID = dplyr::row_number())

## For fitted model
projected_model <- as.matrix(scaled_two_curvy_data_model) %*% projection_scaled

projected_model_df <- projected_model |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::rename(c("proj1" = "...1",
                  "proj2" = "...2")) |>
  dplyr::mutate(ID = dplyr::row_number())

model_df <- dplyr::left_join(
  distance_df_small_edges_two_curvy |> select(-distance),
  projected_model_df,
  by = c("from" = "ID"))

names(model_df)[3:NCOL(model_df)] <- paste0(names(projected_model_df)[-NCOL(projected_model_df)], "_from")

model_df <- dplyr::left_join(model_df, projected_model_df, by = c("to" = "ID"))
names(model_df)[(2 + NCOL(projected_model_df)):NCOL(model_df)] <- paste0(names(projected_model_df)[-NCOL(projected_model_df)], "_to")

## Projection with Delauany triangulation

model_df1 <- dplyr::left_join(
  distance_two_curvy |> select(-distance),
  projected_model_df,
  by = c("from" = "ID"))

names(model_df1)[3:NCOL(model_df1)] <- paste0(names(projected_model_df)[-NCOL(projected_model_df)], "_from")

model_df1 <- dplyr::left_join(model_df1, projected_model_df, by = c("to" = "ID"))
names(model_df1)[(2 + NCOL(projected_model_df)):NCOL(model_df1)] <- paste0(names(projected_model_df)[-NCOL(projected_model_df)], "_to")

## With benchmark option 2

model_df2 <- dplyr::left_join(
  distance_df_small_edges_two_curvy3 |> select(-distance),
  projected_model_df,
  by = c("from" = "ID"))

names(model_df2)[3:NCOL(model_df2)] <- paste0(names(projected_model_df)[-NCOL(projected_model_df)], "_from")

model_df2 <- dplyr::left_join(model_df2, projected_model_df, by = c("to" = "ID"))
names(model_df2)[(2 + NCOL(projected_model_df)):NCOL(model_df2)] <- paste0(names(projected_model_df)[-NCOL(projected_model_df)], "_to")

axes_obj <- gen_axes(
  proj = projection * 5,
  limits = 0.5,
  axis_pos_x = -0.4,
  axis_pos_y = -0.4,
  axis_labels = names(scaled_two_curvy_data),
  threshold = 0.03)

axes <- axes_obj$axes
circle <- axes_obj$circle

two_curvy_proj_first_model1_umap <- projected_df |>
  ggplot(
    aes(
      x = proj1,
      y = proj2)) +
  geom_segment(
    data = model_df,
    aes(
      x = proj1_from,
      y = proj2_from,
      xend = proj1_to,
      yend = proj2_to),
    color = "#000000") + #31a354
  geom_point(
    size = 0.5,
    alpha = 0.2,
    color = "#636363") +
  geom_segment(
    data=axes,
    aes(x=x1, y=y1, xend=x2, yend=y2),
    colour="grey70") +
  geom_text(
    data=axes,
    aes(x=x2, y=y2),
    label=rownames(axes),
    colour="grey50",
    size = 4) +
  geom_path(
    data=circle,
    aes(x=c1, y=c2), colour="grey70") +
  coord_fixed() +
  xlim(c(-0.5, 0.5)) +
  ylim(c(-0.5, 0.5)) +
  interior_annotation("c2", sc_ltr_pos) +
  theme(
    legend.position = "none"
  )

## Projection with Delauany triangulation

two_curvy_proj_model_delaunay_umap <- projected_df |>
  ggplot(
    aes(
      x = proj1,
      y = proj2)) +
  geom_segment(
    data = model_df1,
    aes(
      x = proj1_from,
      y = proj2_from,
      xend = proj1_to,
      yend = proj2_to),
    color = "#000000") + #31a354
  geom_point(
    size = 0.5,
    alpha = 0.2,
    color = "#636363") +
  geom_segment(
    data=axes,
    aes(x=x1, y=y1, xend=x2, yend=y2),
    colour="grey70") +
  geom_text(
    data=axes,
    aes(x=x2, y=y2),
    label=rownames(axes),
    colour="grey50",
    size = 4) +
  geom_path(
    data=circle,
    aes(x=c1, y=c2), colour="grey70") +
  coord_fixed() +
  xlim(c(-0.5, 0.5)) +
  ylim(c(-0.5, 0.5)) +
  interior_annotation("a2", sc_ltr_pos) +
  theme(
    legend.position = "none"
  )

two_curvy_proj_model_benchmark2_umap <- projected_df |>
  ggplot(
    aes(
      x = proj1,
      y = proj2)) +
  geom_segment(
    data = model_df2,
    aes(
      x = proj1_from,
      y = proj2_from,
      xend = proj1_to,
      yend = proj2_to),
    color = "#000000") + #31a354
  geom_point(
    size = 0.5,
    alpha = 0.2,
    color = "#636363") +
  geom_segment(
    data=axes,
    aes(x=x1, y=y1, xend=x2, yend=y2),
    colour="grey70") +
  geom_text(
    data=axes,
    aes(x=x2, y=y2),
    label=rownames(axes),
    colour="grey50",
    size = 4) +
  geom_path(
    data=circle,
    aes(x=c1, y=c2), colour="grey70") +
  coord_fixed() +
  xlim(c(-0.5, 0.5)) +
  ylim(c(-0.5, 0.5)) +
  interior_annotation("b2", sc_ltr_pos) +
  theme(
    legend.position = "none"
  )

gen_tuning_lg_plots <- function() {

  #wrap_elements(
  free(distance_points_umap) + # + ggtitle("(a)") + theme(plot.title.position = "plot", plot.title = element_text(hjust = -0.3))) /
    wrap_plots(trimesh_two_curvy_umap, trimesh_two_curvy_removed2_umap,
               trimesh_two_curvy_removed1_umap, two_curvy_proj_model_delaunay_umap,
               two_curvy_proj_model_benchmark2_umap, two_curvy_proj_first_model1_umap,
               ncol = 3) +
    plot_layout(heights = c(0.8, 1))

}
