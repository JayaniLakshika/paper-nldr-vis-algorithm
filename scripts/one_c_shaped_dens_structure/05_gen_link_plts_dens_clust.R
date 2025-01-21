### This script is to generate error plot for one_c_shaped_dens and one_c_shaped_uni_dens cluster with tSNE
library(readr)
library(dplyr)
library(quollr)
library(ggplot2)
library(colorspace)
library(patchwork)
set.seed(20240110)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)

# creating Standardization function
standardize = function(x){
  z <- (x - mean(x)) / sd(x)
  return( z)
}

interior_annotation <- function(label, position = c(0.92, 0.92), cex = 1, col="grey70") {
  annotation_custom(grid::textGrob(label = label,
                                   x = unit(position[1], "npc"), y = unit(position[2], "npc"),
                                   gp = grid::gpar(cex = cex, col=col)))
}

##1. With one_c_shaped_dens

one_c_shaped_data <- read_rds(here::here("data/one_c_shaped_dens_structure/one_c_shaped_dens_data.rds"))
one_c_shaped_data <- one_c_shaped_data |>
  mutate(ID = row_number())

tsne_one_c_shaped <- read_rds(file = "data/one_c_shaped_dens_structure/one_c_shaped_dens_structure_tsne_perplexity_52.rds")

tsne_one_c_shaped_scaled_obj <- gen_scaled_data(
  data = tsne_one_c_shaped)

tsne_one_c_shaped_scaled <- tsne_one_c_shaped_scaled_obj$scaled_nldr |>
  mutate(ID = row_number())
lim1 <- tsne_one_c_shaped_scaled_obj$lim1
lim2 <- tsne_one_c_shaped_scaled_obj$lim2
r2 <- diff(lim2)/diff(lim1)

num_bins <- 20

## hexagon binning to have regular hexagons
hb_obj_one_c_shaped <- hex_binning(
  data = tsne_one_c_shaped_scaled,
  bin1 = num_bins,
  r2 = r2)

a1_2 <- calc_bins_y(
  bin1 = num_bins,
  r2 = r2
)$a1

## Data set with all centroids
all_centroids_df2 <- hb_obj_one_c_shaped$centroids

## Generate all coordinates of hexagons
hex_grid2 <- hb_obj_one_c_shaped$hex_poly

## To obtain the standardise counts within hexbins
counts_df2 <- hb_obj_one_c_shaped$std_cts
df_bin_centroids2 <- extract_hexbin_centroids(
  centroids_df = all_centroids_df2,
  counts_df = counts_df2) |>
  filter(drop_empty == FALSE)

tsne_data_one_c_shaped_with_hb_id <- hb_obj_one_c_shaped$data_hb_id
df_all_one_curvy2 <- dplyr::bind_cols(one_c_shaped_data |> select(-ID),
                                      tsne_data_one_c_shaped_with_hb_id)

df_bin_one_curvy2 <- avg_highd_data(data = df_all_one_curvy2, col_start = "x")

## Compute error
error_df_one_curvy_abs <- augment(
  df_bin_centroids = df_bin_centroids2,
  df_bin = df_bin_one_curvy2,
  training_data = one_c_shaped_data,
  newdata = NULL,
  type_NLDR = "tSNE",
  col_start = "x")

error_df_one_curvy_abs <- error_df_one_curvy_abs |>
  bind_cols(tsne_one_c_shaped_scaled |>
              select(-ID))

error_df_one_curvy_abs <- error_df_one_curvy_abs |>
  mutate(sqrt_row_wise_total_error = sqrt(row_wise_total_error)) |>
  mutate(sqrt_row_wise_total_error = standardize(sqrt_row_wise_total_error))

quant_val <- quantile(error_df_one_curvy_abs$sqrt_row_wise_total_error,
                      probs = seq(0, 1, 0.12))

error_df_one_curvy_abs <- error_df_one_curvy_abs |>
  mutate(error_cat = if_else(
    sqrt_row_wise_total_error <= quant_val[1], "first", if_else(
      sqrt_row_wise_total_error <= quant_val[2], "second", if_else(
        sqrt_row_wise_total_error <= quant_val[3], "third",
        if_else(
          sqrt_row_wise_total_error <= quant_val[4], "fourth",
          if_else(
            sqrt_row_wise_total_error <= quant_val[5], "fifth",
            if_else(
              sqrt_row_wise_total_error <= quant_val[6], "sixth", if_else(
                sqrt_row_wise_total_error <= quant_val[7], "seventh",
                if_else(
                  sqrt_row_wise_total_error <= quant_val[8], "eighth",
                  "nineth")))))))))

error_plot_tsne <- error_df_one_curvy_abs |>
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             colour = factor(error_cat,
                             levels = c("first", "second", "third",
                                        "fourth", "fifth", "sixth",
                                        "seventh", "eighth", "nineth")))) +
  geom_point(alpha=0.7) +
  scale_color_manual(values=c('#ffffcc','#ffeda0','#fed976','#feb24c',
                              '#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')) +
  theme(
    aspect.ratio = 1
  ) +
  interior_annotation("b1",
                      position = c(0.08, 0.95),
                      cex = 1.5)


## Triangulate bin centroids
tr1_object_c_shaped_structure <- tri_bin_centroids(
  df_bin_centroids2, x = "c_x", y = "c_y")
tr_from_to_df_c_shaped_structure <- gen_edges(
  tri_object = tr1_object_c_shaped_structure)

## Compute 2D distances
distance_c_shaped_structure <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_c_shaped_structure,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_c_shaped_structure <- find_lg_benchmark(
  distance_edges = distance_c_shaped_structure,
  distance_col = "distance")

tr_df <- distinct(tibble::tibble(
  x = c(tr_from_to_df_c_shaped_structure[["x_from"]], tr_from_to_df_c_shaped_structure[["x_to"]]),
  y = c(tr_from_to_df_c_shaped_structure[["y_from"]], tr_from_to_df_c_shaped_structure[["y_to"]])))

distance_df_small_edges_c_shaped_structure <- distance_c_shaped_structure |>
  filter(distance < benchmark_c_shaped_structure)

tr_from_to_df_c_shaped_structure <- inner_join(
  tr_from_to_df_c_shaped_structure, distance_df_small_edges_c_shaped_structure,
  by = c("from", "to"))

trimesh_removed_c_shaped_structure <- ggplot() +
  geom_segment(data = tr_from_to_df_c_shaped_structure,
               aes(
                 x = x_from,
                 y = y_from,
                 xend = x_to,
                 yend = y_to),
               colour = "#33a02c",
               linewidth = 1) +
  geom_point(data = tsne_one_c_shaped_scaled,
             aes(
               x = tSNE1,
               y = tSNE2
             ),
             alpha=0.1) +
  theme(aspect.ratio = 1)


### Define type column
df <- df_all_one_curvy2 |>
  dplyr::select(tidyselect::starts_with("x")) |>
  dplyr::mutate(type = error_df_one_curvy_abs$error_cat) ## original dataset

df_b <- df_bin_one_curvy2 |>
  dplyr::filter(hb_id %in% df_bin_centroids2$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids2$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

df_exe <- dplyr::bind_rows(df_b, df)

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges_c_shaped_structure$from,
                         lineTo = distance_df_small_edges_c_shaped_structure$to,
                         group = factor(df_exe$type,
                                        c("first", "second", "third",
                                          "fourth", "fifth", "sixth",
                                          "seventh", "eighth", "nineth","model")), pointSize = append(rep(1, NROW(df_b)), rep(2, NROW(df))),
                         levelColors = c('#ffffcc','#ffeda0','#fed976','#feb24c',
                                         '#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026',
                                         "#000000"))

df_all_one_curvy2 <- df_all_one_curvy2 |> dplyr::select(-ID, -hb_id)

show_link_plots(df_all = df_all_one_curvy2, df_b = df_bin_one_curvy2,
                df_b_with_center_data = df_bin_centroids2, benchmark_value = benchmark_c_shaped_structure,
                distance = distance_c_shaped_structure, distance_col = "distance",
                use_default_benchmark_val = FALSE, col_start = "x", type_nldr = "tSNE", r2 = r2)

