### This script is to generate error three bin choices
library(readr)
library(quollr)
library(dplyr)

data_two_curvy <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_data.rds")
data_two_curvy <- data_two_curvy |>
  mutate(ID = row_number())

tsne_two_curvy <- read_rds(file = "data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_47.rds") |>
  mutate(ID = row_number())


### Option 1

two_curvy_model1 <- fit_highd_model(
  highd_data = data_two_curvy,
  nldr_data = tsne_two_curvy,
  bin1 = 15,
  q = 0.1,
  benchmark_highdens = 0)

df_bin_centroids_two_curvy1 <- two_curvy_model1$model_2d
df_bin_two_curvy1 <- two_curvy_model1$model_highd

error_rm_two_curvy1 <- data.frame(matrix(nrow = 0, ncol = 0))

## To initialize benchmark values to remove low density hexagons
benchmark_rm_hex_vec <- seq(0, 0.47, by=0.03)

for (benchmark_rm_lwd in benchmark_rm_hex_vec) {

  df_bin_centroids_two_curvy_high_dens <- df_bin_centroids_two_curvy1 |>
    filter(std_counts > benchmark_rm_lwd)

  df_bin_two_curvy_high_dens <- df_bin_two_curvy1 |>
    filter(hexID %in% df_bin_centroids_two_curvy_high_dens$hexID)

  ## Compute error
  error_df <- glance(
    model_2d = df_bin_centroids_two_curvy_high_dens,
    model_highd = df_bin_two_curvy_high_dens,
    highd_data = data_two_curvy) |>
    mutate(benchmark_rm_lwd = round(benchmark_rm_lwd, 2),
           a1 = paste0("a[1] == ", round(two_curvy_model1$hb_obj$a1, 2)),
           bin1 = two_curvy_model1$hb_obj$bins[1],
           bin2 = two_curvy_model1$hb_obj$bins[2],
           b = bin1 * bin2,
           m = NROW(df_bin_centroids_two_curvy_high_dens),
           mean_counts = sum(df_bin_centroids_two_curvy_high_dens$bin_counts)/NROW(df_bin_centroids_two_curvy_high_dens))

  error_rm_two_curvy1 <- bind_rows(error_rm_two_curvy1, error_df)

}

### Option 2

two_curvy_model2 <- fit_highd_model(
  highd_data = data_two_curvy,
  nldr_data = tsne_two_curvy,
  bin1 = 24,
  q = 0.1,
  benchmark_highdens = 0)

df_bin_centroids_two_curvy2 <- two_curvy_model2$model_2d
df_bin_two_curvy2 <- two_curvy_model2$model_highd

error_rm_two_curvy2 <- data.frame(matrix(nrow = 0, ncol = 0))

## To initialize benchmark values to remove low density hexagons
benchmark_rm_hex_vec <- seq(0, 0.47, by=0.03)

for (benchmark_rm_lwd in benchmark_rm_hex_vec) {

  df_bin_centroids_two_curvy_high_dens <- df_bin_centroids_two_curvy2 |>
    filter(std_counts > benchmark_rm_lwd)

  df_bin_two_curvy_high_dens <- df_bin_two_curvy2 |>
    filter(hexID %in% df_bin_centroids_two_curvy_high_dens$hexID)

  ## Compute error
  error_df <- glance(
    model_2d = df_bin_centroids_two_curvy_high_dens,
    model_highd = df_bin_two_curvy_high_dens,
    highd_data = data_two_curvy) |>
    mutate(benchmark_rm_lwd = round(benchmark_rm_lwd, 2),
           a1 = paste0("a[1] == ", round(two_curvy_model2$hb_obj$a1, 2)),
           bin1 = two_curvy_model2$hb_obj$bins[1],
           bin2 = two_curvy_model2$hb_obj$bins[2],
           b = bin1 * bin2,
           m = NROW(df_bin_centroids_two_curvy_high_dens),
           mean_counts = sum(df_bin_centroids_two_curvy_high_dens$bin_counts)/NROW(df_bin_centroids_two_curvy_high_dens))

  error_rm_two_curvy2 <- bind_rows(error_rm_two_curvy2, error_df)

}

### Option 3

two_curvy_model3 <- fit_highd_model(
  highd_data = data_two_curvy,
  nldr_data = tsne_two_curvy,
  bin1 = 48,
  q = 0.1,
  benchmark_highdens = 0)

df_bin_centroids_two_curvy3 <- two_curvy_model3$model_2d
df_bin_two_curvy3 <- two_curvy_model3$model_highd

error_rm_two_curvy3 <- data.frame(matrix(nrow = 0, ncol = 0))

## To initialize benchmark values to remove low density hexagons
benchmark_rm_hex_vec <- seq(0, 0.47, by=0.03)

for (benchmark_rm_lwd in benchmark_rm_hex_vec) {

  df_bin_centroids_two_curvy_high_dens <- df_bin_centroids_two_curvy3 |>
    filter(std_counts > benchmark_rm_lwd)

  df_bin_two_curvy_high_dens <- df_bin_two_curvy3 |>
    filter(hexID %in% df_bin_centroids_two_curvy_high_dens$hexID)

  ## Compute error
  error_df <- glance(
    model_2d = df_bin_centroids_two_curvy_high_dens,
    model_highd = df_bin_two_curvy_high_dens,
    highd_data = data_two_curvy) |>
    mutate(benchmark_rm_lwd = round(benchmark_rm_lwd, 2),
           a1 = paste0("a[1] == ", round(two_curvy_model3$hb_obj$a1, 2)),
           bin1 = two_curvy_model3$hb_obj$bins[1],
           bin2 = two_curvy_model3$hb_obj$bins[2],
           b = bin1 * bin2,
           m = NROW(df_bin_centroids_two_curvy_high_dens),
           mean_counts = sum(df_bin_centroids_two_curvy_high_dens$bin_counts)/NROW(df_bin_centroids_two_curvy_high_dens))

  error_rm_two_curvy3 <- bind_rows(error_rm_two_curvy3, error_df)

}

error_rm_two_curvy <- bind_rows(
  error_rm_two_curvy1,
  error_rm_two_curvy2,
  error_rm_two_curvy3
)

write_rds(error_rm_two_curvy, here::here("data/two_nonlinear/error_rm_lwd_diff_bin.rds"))
