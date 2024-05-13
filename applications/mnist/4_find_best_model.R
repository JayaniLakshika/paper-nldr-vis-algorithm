library(quollr)
library(dplyr)
library(readr)

## Import data
training_data_mnist <- read_rds("data/mnist/mnist_10_pcs_of_digit_1.rds")
training_data_mnist <- training_data_mnist |>
  mutate(ID = 1:NROW(training_data_mnist))

pacmap_minst <- read_rds("data/mnist/mnist_pacmap.rds")
mnist_scaled_obj <- gen_scaled_data(
  data = pacmap_minst)
pacmap_minst_scaled <- mnist_scaled_obj$scaled_nldr

## Compute range of upper bound
lim1 <- mnist_scaled_obj$lim1
lim2 <- mnist_scaled_obj$lim2
r2_mnist <- diff(lim2)/diff(lim1)

## To initialise number of bins along the x-axis
bin1_vec <- 2:50

error_minst <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec) {

 bin2 <- calc_bins_y(bin1 = xbins, q = 0.1, r2 = r2_mnist)$bin2

  mnist_model <- fit_highd_model(
    training_data = training_data_mnist,
    emb_df = pacmap_minst_scaled,
    bin1 = xbins,
    q = 0.1,
    r2 = r2_mnist,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "PC"
  )

  df_bin_centroids_mnist <- mnist_model$df_bin_centroids
  df_bin_mnist <- mnist_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_mnist,
    df_bin = df_bin_mnist,
    training_data = training_data_mnist,
    newdata = NULL,
    type_NLDR = "PaCMAP",
    col_start = "PC") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_mnist))

  error_minst <- bind_rows(error_minst, error_df)

}

error_minst |>
  arrange(MSE)


ggplot(error_minst, aes(x = b_non_empty,
                     y = log(MSE))) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 341, linetype="solid",
             color = "black", linewidth=0.8, alpha = 0.5)
