library(quollr)
library(dplyr)
library(readr)
library(ggplot2)

## Import data
training_data_gau <- read_rds("data/five_gau_clusters/data_five_gau_training.rds")

tsne_gau <- read_rds("data/five_gau_clusters/tsne_data_five_gau_61.rds")
gau_scaled_obj <- gen_scaled_data(
  data = tsne_gau)
tsne_gau_scaled <- gau_scaled_obj$scaled_nldr

## To initialise number of bins along the x-axis
bin1_vec <- 2:61

lim1 <- gau_scaled_obj$lim1
lim2 <- gau_scaled_obj$lim2
r2_gau <- diff(lim2)/diff(lim1)

error_gau <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec) {

  bin2 <- calc_bins_y(bin1 = xbins, r2 = r2_gau)$bin2

  gau_model <- fit_highd_model(
    training_data = training_data_gau,
    emb_df = tsne_gau_scaled,
    bin1 = xbins,
    r2 = r2_gau,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "x"
  )

  df_bin_centroids_gau <- gau_model$df_bin_centroids
  df_bin_gau <- gau_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_gau,
    df_bin = df_bin_gau,
    training_data = training_data_gau,
    newdata = NULL,
    type_NLDR = "tSNE",
    col_start = "x") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_gau))

  error_gau <- bind_rows(error_gau, error_df)

}

error_gau |>
  arrange(MSE)


ggplot(error_gau, aes(x = b_non_empty,
                      y = log(MSE))) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 65, linetype="solid",
             color = "black", linewidth=0.8, alpha = 0.5)
