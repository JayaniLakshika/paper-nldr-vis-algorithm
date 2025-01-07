### This script is to generate error plot for one_c_shaped_dens and one_c_shaped_uni_dens cluster with tSNE
library(readr)
library(dplyr)
library(quollr)
library(ggplot2)
library(colorspace)
library(patchwork)
set.seed(20240110)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

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

one_c_shaped_data <- read_rds(here::here("data/one_c_shaped_dens_structure/one_c_shaped_dens_data.rds"))
one_c_shaped_data <- one_c_shaped_data |>
  mutate(ID = row_number())

##1. With one_c_shaped_dens

tsne_one_c_shaped <- read_rds(file = "data/one_c_shaped_dens_structure/one_c_shaped_dens_structure_tsne_perplexity_52.rds")

tsne_one_c_shaped_scaled_obj <- gen_scaled_data(
  data = tsne_one_c_shaped)

tsne_one_c_shaped_scaled <- tsne_one_c_shaped_scaled_obj$scaled_nldr |>
  mutate(ID = row_number())
lim1 <- tsne_one_c_shaped_scaled_obj$lim1
lim2 <- tsne_one_c_shaped_scaled_obj$lim2
r2 <- diff(lim2)/diff(lim1)

num_bins <- 17

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
df_all_one_curvy2 <- dplyr::bind_cols(one_c_shaped_data,
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
  mutate(row_wise_total_error = standardize(row_wise_total_error))

error_plot_tsne <- error_df_one_curvy_abs |>
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             colour = row_wise_total_error)) +
  geom_point(alpha=0.5) +
  scale_colour_continuous_sequential(palette = "YlOrRd", n_interp = 20) +
  theme(
    aspect.ratio = 1
  ) +
  interior_annotation("a1",
                      position = c(0.08, 0.95),
                      cex = 1.5)

##2. With one_c_shaped_uni_dens

one_c_shaped_data <- read_rds(here::here("data/one_c_shaped_dens_structure/one_c_shaped_uni_dens_data.rds"))
one_c_shaped_data <- one_c_shaped_data |>
  mutate(ID = row_number())

tsne_one_c_shaped <- read_rds(file = "data/one_c_shaped_dens_structure/one_c_shaped_uni_dens_structure_tsne_perplexity_52.rds")

tsne_one_c_shaped_scaled_obj <- gen_scaled_data(
  data = tsne_one_c_shaped)

tsne_one_c_shaped_scaled <- tsne_one_c_shaped_scaled_obj$scaled_nldr |>
  mutate(ID = row_number())
lim1 <- tsne_one_c_shaped_scaled_obj$lim1
lim2 <- tsne_one_c_shaped_scaled_obj$lim2
r2 <- diff(lim2)/diff(lim1)

num_bins <- 13

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
df_all_one_curvy2 <- dplyr::bind_cols(one_c_shaped_data,
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
  mutate(row_wise_total_error = standardize(row_wise_total_error))

error_plot_tsne_uni <- error_df_one_curvy_abs |>
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             colour = row_wise_total_error)) +
  geom_point(alpha=0.5) +
  scale_colour_continuous_sequential(palette = "YlOrRd", n_interp = 20) +
  theme(
    aspect.ratio = 1
  ) +
  interior_annotation("b1",
                      position = c(0.08, 0.95),
                      cex = 1.5)



## Function to print all the plots together

generate_error_plots_one_c_shaped <- function(){

  error_plot_tsne + error_plot_tsne_uni +
    plot_layout(guides = "collect", ncol = 2) &
    theme(legend.position='none')

}

