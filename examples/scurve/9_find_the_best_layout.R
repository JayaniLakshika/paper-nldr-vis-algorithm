## Import necessary libraries
library(quollr)
library(dplyr)
library(reader)

## Import S-curve data
training_data_scurve <- read_rds("data/s_curve/s_curve_training.rds")

## Import layouts data
tsne_scurve <- read_rds(file = "data/s_curve/s_curve_tsne_27.rds") |>
  mutate(method = "tSNE") |>
  rename(c("emb1" = "tSNE1",
           "emb2" = "tSNE2"))

umap_scurve <- read_rds(file = "data/s_curve/s_curve_umap.rds") |>
  mutate(method = "UMAP") |>
  rename(c("emb1" = "UMAP1",
           "emb2" = "UMAP2"))

phate_scurve <- read_rds(file = "data/s_curve/s_curve_phate.rds") |>
  mutate(method = "PHATE") |>
  rename(c("emb1" = "PHATE1",
           "emb2" = "PHATE2"))

trimap_scurve <- read_rds(file = "data/s_curve/s_curve_trimap.rds") |>
  mutate(method = "TriMAP") |>
  rename(c("emb1" = "TriMAP1",
           "emb2" = "TriMAP2"))

pacmap_scurve <- read_rds(file = "data/s_curve/s_curve_pacmap.rds") |>
  mutate(method = "PaCMAP") |>
  rename(c("emb1" = "PaCMAP1",
           "emb2" = "PaCMAP2"))

nldr_scurve <- bind_rows(tsne_scurve, umap_scurve, phate_scurve,
                         trimap_scurve, pacmap_scurve)

## To initialize number of bins along the x-axis
bin1_vec_scurve <- 2:19 #sqrt(NROW(training_data_scurve)/2)

error_scurve <- data.frame(matrix(nrow = 0, ncol = 0))

for (nldr in c("tSNE", "UMAP", "PHATE", "PaCMAP", "TriMAP")) {

  nldr_scurve_filtered <- nldr_scurve |>
    filter(method == nldr)

  scurve_scaled_obj <- gen_scaled_data(
    data = nldr_scurve_filtered)

  nldr_scurve_scaled <- scurve_scaled_obj$scaled_nldr
  lim1 <- scurve_scaled_obj$lim1
  lim2 <- scurve_scaled_obj$lim2
  r2 <- diff(lim2)/diff(lim1)

  for (xbins in bin1_vec_scurve) {

    bin2 <- calc_bins_y(bin1 = xbins, r2 = r2)$bin2

    scurve_model <- fit_highd_model(
      training_data = training_data_scurve,
      emb_df = nldr_scurve_scaled,
      bin1 = xbins,
      r2 = r2,
      is_bin_centroid = TRUE,
      is_rm_lwd_hex = FALSE,
      col_start_highd = "x"
    )

    df_bin_centroids_scurve <- scurve_model$df_bin_centroids
    df_bin_scurve <- scurve_model$df_bin

    ## Compute error
    error_df <- glance(
      df_bin_centroids = df_bin_centroids_scurve,
      df_bin = df_bin_scurve,
      training_data = training_data_scurve,
      newdata = NULL,
      type_NLDR = nldr,
      col_start = "x") |>
      mutate(bin1 = xbins,
             bin2 = bin2,
             b = bin1 * bin2,
             b_non_empty = NROW(df_bin_centroids_scurve),
             method = nldr)

    error_scurve <- bind_rows(error_scurve, error_df)

  }


}


ggplot(error_scurve,
       aes(x = b_non_empty,
           y = log(Error),
           group = method,
           colour = method)) +
  geom_point(size = 1) +
  geom_line() +
  # geom_vline(xintercept = 190, linetype="solid",
  #            color = "black", linewidth=0.8, alpha = 0.5) +
  ylab("log(MSE)") +
  xlab("b") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10, angle = 90))
## Selected PaCMAP with b_non_empty = 42 (11 * 12 = 132)
