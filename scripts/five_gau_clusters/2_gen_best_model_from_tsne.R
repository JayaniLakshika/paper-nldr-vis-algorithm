library(quollr)
library(dplyr)
library(readr)
library(ggplot2)

## Import data
training_data_gau <- read_rds("data/five_gau_clusters/data_five_gau.rds")

tsne_gau <- read_rds("data/five_gau_clusters/tsne_data_five_gau_61.rds")
names(pacmap_gau) <- c("emb1", "emb2", "ID")

gau_scaled_obj <- gen_scaled_data(
  data = tsne_gau)
tsne_gau_scaled <- gau_scaled_obj$scaled_nldr

## To initialise number of bins along the x-axis
bin1_vec <- 2:71 # sqrt(NROW(training_data_gau)/r2_gau)
#buff_vec <- seq(0.05, 0.2, by = 0.01)

lim1 <- gau_scaled_obj$lim1
lim2 <- gau_scaled_obj$lim2
r2_gau <- diff(lim2)/diff(lim1)

error_gau <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec) {
  #for(q in buff_vec) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_gau, q = 0.1)
  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  gau_model <- fit_highd_model(
    highd_data = training_data_gau,
    nldr_data = pacmap_gau_scaled,
    bin1 = xbins,
    r2 = r2_gau,
    is_bin_centroid = TRUE,
    q = 0.1
  )

  df_bin_centroids_gau <- gau_model$df_bin_centroids
  df_bin_gau <- gau_model$df_bin

  ## Compute error
  error_df <- glance(
    model_2d = df_bin_centroids_gau,
    model_highd = df_bin_gau,
    highd_data = training_data_gau) |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_gau),
           a1 = a1)

  error_gau <- bind_rows(error_gau, error_df)

  #}

}

error_gau |>
  arrange(MSE)


ggplot(error_gau, aes(x = a1,
                      y = MSE)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 65, linetype="solid",
             color = "black", linewidth=0.8, alpha = 0.5)
