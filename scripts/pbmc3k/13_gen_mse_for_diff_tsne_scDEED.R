## This script is to generate MSE for hyper-parameters in tSNE suggested by Chen and scDEED

library(readr)
library(quollr)
library(dplyr)

training_data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50_scdeed.rds")
training_data_pbmc <- training_data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(training_data_pbmc))

names(training_data_pbmc) <- append(paste0("x", 1:9), "ID")

## For tsne
tsne_pbmc <- read_rds("data/pbmc3k/pbmc_scdeed_tsne_perplexity_30.rds")
tsne_pbmc <- as_tibble(tsne_pbmc)
names(tsne_pbmc) <- c("emb1", "emb2")
tsne_pbmc <- tsne_pbmc |>
  mutate(ID = 1:NROW(tsne_pbmc))

pbmc_scaled_obj_tsne <- gen_scaled_data(
  data = tsne_pbmc)
tsne_pbmc_scaled <- pbmc_scaled_obj_tsne$scaled_nldr

lim1 <- pbmc_scaled_obj_tsne$lim1
lim2 <- pbmc_scaled_obj_tsne$lim2
r2_tsne <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:78 #sqrt(NROW(training_data_pbmc)/r2_tsne)

error_pbmc_tsne <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_tsne, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    highd_data = training_data_pbmc,
    nldr_data = tsne_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_tsne,
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
           method = "tsne_perplexity_30",
           a1 = a1)

  error_pbmc_tsne <- bind_rows(error_pbmc_tsne, error_df)

}

write_rds(error_pbmc_tsne, "data/pbmc3k/error_scdeed_pbmc_tsne_perplexity_30.rds")

#################################

## For tsne
tsne_pbmc <- read_rds("data/pbmc3k/pbmc_scdeed_tsne_perplexity_320.rds")
tsne_pbmc <- as_tibble(tsne_pbmc)
names(tsne_pbmc) <- c("emb1", "emb2")
tsne_pbmc <- tsne_pbmc |>
  mutate(ID = 1:NROW(tsne_pbmc))

pbmc_scaled_obj_tsne <- gen_scaled_data(
  data = tsne_pbmc)
tsne_pbmc_scaled <- pbmc_scaled_obj_tsne$scaled_nldr

lim1 <- pbmc_scaled_obj_tsne$lim1
lim2 <- pbmc_scaled_obj_tsne$lim2
r2_tsne <- diff(lim2)/diff(lim1)

bin1_vec_pbmc <- 2:87 #sqrt(NROW(training_data_pbmc)/r2_tsne)

error_pbmc_tsne <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_tsne, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    highd_data = training_data_pbmc,
    nldr_data = tsne_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_tsne,
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
           method = "tsne_perplexity_320",
           a1 = a1)

  error_pbmc_tsne <- bind_rows(error_pbmc_tsne, error_df)

}

write_rds(error_pbmc_tsne, "data/pbmc3k/error_scdeed_pbmc_tsne_perplexity_320.rds")

