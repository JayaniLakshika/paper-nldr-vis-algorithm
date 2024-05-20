library(quollr)
library(dplyr)
library(readr)
library(ggplot2)

## Import data
training_data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50.rds")
training_data_pbmc <- training_data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(training_data_pbmc))

umap_pbmc <- read_rds("data/pbmc3k/pbmc_umap_30_min_dist_0.3.rds")
pbmc_scaled_obj <- gen_scaled_data(
  data = umap_pbmc)
umap_pbmc_scaled <- pbmc_scaled_obj$scaled_nldr

## To initialise number of bins along the x-axis
bin1_vec <- 2:50

lim1 <- pbmc_scaled_obj$lim1
lim2 <- pbmc_scaled_obj$lim2
r2_pbmc <- diff(lim2)/diff(lim1)

error_pbmc <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec) {

  bin2 <- calc_bins_y(bin1 = xbins, r2 = r2_pbmc)$bin2

  pbmc_model <- fit_highd_model(
    training_data = training_data_pbmc,
    emb_df = umap_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_pbmc,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "PC_"
  )

  df_bin_centroids_pbmc <- pbmc_model$df_bin_centroids
  df_bin_pbmc <- pbmc_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_pbmc,
    df_bin = df_bin_pbmc,
    training_data = training_data_pbmc,
    newdata = NULL,
    type_NLDR = "UMAP",
    col_start = "PC_") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_pbmc))

  error_pbmc <- bind_rows(error_pbmc, error_df)

}

error_pbmc |>
  arrange(MSE)


ggplot(error_pbmc, aes(x = b_non_empty,
                        y = log(MSE))) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 114, linetype="solid",
             color = "black", linewidth=0.8, alpha = 0.5)
