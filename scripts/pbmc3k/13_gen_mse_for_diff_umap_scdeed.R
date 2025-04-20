## This script is to generate MSE for hyper-parameter suggested by Chen and scDEED
library(readr)
library(quollr)
library(dplyr)

training_data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50_scdeed.rds")
training_data_pbmc <- training_data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(training_data_pbmc))

names(training_data_pbmc) <- append(paste0("x", 1:9), "ID")

## For umap
umap_pbmc <- read_rds("data/pbmc3k/pbmc_scdeed_umap_n_neighbors_30_min_dist_0.3.rds")
umap_pbmc <- as_tibble(umap_pbmc)
names(umap_pbmc) <- c("emb1", "emb2")
umap_pbmc <- umap_pbmc |>
  mutate(ID = 1:NROW(umap_pbmc))

pbmc_scaled_obj_umap <- gen_scaled_data(
  data = umap_pbmc)
umap_pbmc_scaled <- pbmc_scaled_obj_umap$scaled_nldr

lim1 <- pbmc_scaled_obj_umap$lim1
lim2 <- pbmc_scaled_obj_umap$lim2
r2_umap <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:84 #sqrt(NROW(training_data_pbmc)/r2_umap)

error_pbmc_umap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    highd_data = training_data_pbmc,
    nldr_data = umap_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_umap,
    q = 0.1,
    is_bin_centroid = TRUE
  )

  df_bin_centroids_pbmc <- pbmc_model$df_bin_centroids
  df_bin_pbmc <- pbmc_model$df_bin

  ## Compute error
  error_df <- glance(
    model_2d = df_bin_centroids_pbmc,
    model_highd = df_bin_pbmc,
    highd_data = training_data_pbmc) |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_pbmc),
           method = "UMAP_30_min_dist_0.3",
           a1 = a1)

  error_pbmc_umap <- bind_rows(error_pbmc_umap, error_df)

}

write_rds(error_pbmc_umap, "data/pbmc3k/error_scdeed_pbmc_umap_30_min_dist_0.3.rds")

#################################

## For umap
umap_pbmc <- read_rds("data/pbmc3k/pbmc_scdeed_umap_n_neighbors_80_min_dist_0.5.rds")
umap_pbmc <- as_tibble(umap_pbmc)
names(umap_pbmc) <- c("emb1", "emb2")
umap_pbmc <- umap_pbmc |>
  mutate(ID = 1:NROW(umap_pbmc))

pbmc_scaled_obj_umap <- gen_scaled_data(
  data = umap_pbmc)
umap_pbmc_scaled <- pbmc_scaled_obj_umap$scaled_nldr

lim1 <- pbmc_scaled_obj_umap$lim1
lim2 <- pbmc_scaled_obj_umap$lim2
r2_umap <- diff(lim2)/diff(lim1)

bin1_vec_pbmc <- 2:69 #sqrt(NROW(training_data_pbmc)/r2_umap)

error_pbmc_umap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    highd_data = training_data_pbmc,
    nldr_data = umap_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_umap,
    q = 0.1,
    is_bin_centroid = TRUE
  )

  df_bin_centroids_pbmc <- pbmc_model$df_bin_centroids
  df_bin_pbmc <- pbmc_model$df_bin

  ## Compute error
  error_df <- glance(
    model_2d = df_bin_centroids_pbmc,
    model_highd = df_bin_pbmc,
    highd_data = training_data_pbmc) |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_pbmc),
           method = "UMAP_80_min_dist_0.5",
           a1 = a1)

  error_pbmc_umap <- bind_rows(error_pbmc_umap, error_df)

}

write_rds(error_pbmc_umap, "data/pbmc3k/error_scdeed_pbmc_umap_80_min_dist_0.5.rds")

