### This script is to generate error three bin choices
library(readr)
library(quollr)
library(dplyr)

training_data_two_curvy <- read_rds("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_data.rds")
training_data_two_curvy <- training_data_two_curvy |>
  mutate(ID = row_number())

tsne_two_curvy <- read_rds(file = "data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_30.rds")

two_curvy_scaled_obj <- gen_scaled_data(
  data = tsne_two_curvy)

tsne_two_curvy_scaled <- two_curvy_scaled_obj$scaled_nldr |>
  mutate(ID = row_number())
lim1 <- two_curvy_scaled_obj$lim1
lim2 <- two_curvy_scaled_obj$lim2
r2 <- diff(lim2)/diff(lim1)


### Option 1

two_curvy_model1 <- fit_highd_model(
  training_data = training_data_two_curvy,
  emb_df = tsne_two_curvy_scaled,
  bin1 = 15,
  r2 = r2,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_two_curvy1 <- two_curvy_model1$df_bin_centroids
df_bin_two_curvy1 <- two_curvy_model1$df_bin

error_rm_two_curvy1 <- data.frame(matrix(nrow = 0, ncol = 0))

## To initialize benchmark values to remove low density hexagons
benchmark_rm_hex_vec <- seq(0, 0.47, by=0.03)

for (benchmark_rm_lwd in benchmark_rm_hex_vec) {

  df_bin_centroids_two_curvy_high_dens <- df_bin_centroids_two_curvy1 |>
    filter(std_counts > benchmark_rm_lwd)

  df_bin_two_curvy_high_dens <- df_bin_two_curvy1 |>
    filter(hb_id %in% df_bin_centroids_two_curvy_high_dens$hexID)

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_two_curvy_high_dens,
    df_bin = df_bin_two_curvy_high_dens,
    training_data = training_data_two_curvy,
    newdata = NULL,
    type_NLDR = "tsne",
    col_start = "x") |>
    mutate(benchmark_rm_lwd = round(benchmark_rm_lwd, 2),
           bin1 = 10,
           bin2 = 8,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_two_curvy_high_dens),
           mean_counts = sum(df_bin_centroids_two_curvy_high_dens$bin_counts)/NROW(df_bin_centroids_two_curvy_high_dens))

  error_rm_two_curvy1 <- bind_rows(error_rm_two_curvy1, error_df)

}

### Option 2

two_curvy_model2 <- fit_highd_model(
  training_data = training_data_two_curvy,
  emb_df = tsne_two_curvy_scaled,
  bin1 = 23,
  r2 = r2,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_two_curvy2 <- two_curvy_model2$df_bin_centroids
df_bin_two_curvy2 <- two_curvy_model2$df_bin

error_rm_two_curvy2 <- data.frame(matrix(nrow = 0, ncol = 0))

## To initialize benchmark values to remove low density hexagons
benchmark_rm_hex_vec <- seq(0, 0.47, by=0.03)

for (benchmark_rm_lwd in benchmark_rm_hex_vec) {

  df_bin_centroids_two_curvy_high_dens <- df_bin_centroids_two_curvy2 |>
    filter(std_counts > benchmark_rm_lwd)

  df_bin_two_curvy_high_dens <- df_bin_two_curvy2 |>
    filter(hb_id %in% df_bin_centroids_two_curvy_high_dens$hexID)

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_two_curvy_high_dens,
    df_bin = df_bin_two_curvy_high_dens,
    training_data = training_data_two_curvy,
    newdata = NULL,
    type_NLDR = "tsne",
    col_start = "x") |>
    mutate(benchmark_rm_lwd = round(benchmark_rm_lwd, 2),
           bin1 = 13,
           bin2 = 10,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_two_curvy_high_dens),
           mean_counts = sum(df_bin_centroids_two_curvy_high_dens$bin_counts)/NROW(df_bin_centroids_two_curvy_high_dens))

  error_rm_two_curvy2 <- bind_rows(error_rm_two_curvy2, error_df)

}

### Option 3

two_curvy_model3 <- fit_highd_model(
  training_data = training_data_two_curvy,
  emb_df = tsne_two_curvy_scaled,
  bin1 = 48,
  r2 = r2,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_two_curvy3 <- two_curvy_model3$df_bin_centroids
df_bin_two_curvy3 <- two_curvy_model3$df_bin

error_rm_two_curvy3 <- data.frame(matrix(nrow = 0, ncol = 0))

## To initialize benchmark values to remove low density hexagons
benchmark_rm_hex_vec <- seq(0, 0.47, by=0.03)

for (benchmark_rm_lwd in benchmark_rm_hex_vec) {

  df_bin_centroids_two_curvy_high_dens <- df_bin_centroids_two_curvy3 |>
    filter(std_counts > benchmark_rm_lwd)

  df_bin_two_curvy_high_dens <- df_bin_two_curvy3 |>
    filter(hb_id %in% df_bin_centroids_two_curvy_high_dens$hexID)

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_two_curvy_high_dens,
    df_bin = df_bin_two_curvy_high_dens,
    training_data = training_data_two_curvy,
    newdata = NULL,
    type_NLDR = "tsne",
    col_start = "x") |>
    mutate(benchmark_rm_lwd = round(benchmark_rm_lwd, 2),
           bin1 = 23,
           bin2 = 17,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_two_curvy_high_dens),
           mean_counts = sum(df_bin_centroids_two_curvy_high_dens$bin_counts)/NROW(df_bin_centroids_two_curvy_high_dens))

  error_rm_two_curvy3 <- bind_rows(error_rm_two_curvy3, error_df)

}

error_rm_two_curvy <- bind_rows(
  error_rm_two_curvy1,
  error_rm_two_curvy2,
  error_rm_two_curvy3
)

write_rds(error_rm_two_curvy, here::here("data/two_non_linear_diff_shaped_close_clusters/error_rm_lwd_diff_bin.rds"))
