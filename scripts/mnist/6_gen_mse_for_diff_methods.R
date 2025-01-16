library(readr)
library(quollr)
library(dplyr)

training_data_mnist <- read_rds("data/mnist/mnist_10_pcs_of_digit_1.rds")
training_data_mnist <- training_data_mnist |>
  mutate(ID = 1:NROW(training_data_mnist))

## To initialize number of bins along the x-axis
bin1_vec_mnist <- 5:63 #sqrt(NROW(training_data_mnist)/2)

## For umap
umap_mnist <- read_rds("data/mnist/mnist_umap.rds")
mnist_scaled_obj_umap <- gen_scaled_data(
  data = umap_mnist)
umap_mnist_scaled <- mnist_scaled_obj_umap$scaled_nldr

lim1 <- mnist_scaled_obj_umap$lim1
lim2 <- mnist_scaled_obj_umap$lim2
r2_umap <- diff(lim2)/diff(lim1)

error_mnist_umap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_mnist) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  mnist_model <- fit_highd_model(
    training_data = training_data_mnist,
    emb_df = umap_mnist_scaled,
    bin1 = xbins,
    r2 = r2_umap,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = TRUE,
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
    type_NLDR = "UMAP",
    col_start = "PC") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_mnist),
           method = "UMAP",
           a1 = a1)

  error_mnist_umap <- bind_rows(error_mnist_umap, error_df)

}

write_rds(error_mnist_umap, "data/mnist/error_mnist_umap.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_mnist <- 5:63 #sqrt(NROW(training_data_mnist)/2)

## For tsne
#tsne_mnist <- read_rds("data/mnist/mnist_tsne89.rds")
tsne_mnist <- read_rds("data/mnist/mnist_tsne30.rds")
mnist_scaled_obj_tsne <- gen_scaled_data(
  data = tsne_mnist)
tsne_mnist_scaled <- mnist_scaled_obj_tsne$scaled_nldr

lim1 <- mnist_scaled_obj_tsne$lim1
lim2 <- mnist_scaled_obj_tsne$lim2
r2_tsne <- diff(lim2)/diff(lim1)

error_mnist_tsne <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_mnist) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_tsne, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  mnist_model <- fit_highd_model(
    training_data = training_data_mnist,
    emb_df = tsne_mnist_scaled,
    bin1 = xbins,
    r2 = r2_tsne,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = TRUE,
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
    type_NLDR = "tSNE",
    col_start = "PC") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_mnist),
           method = "tSNE",
           a1 = a1)

  error_mnist_tsne <- bind_rows(error_mnist_tsne, error_df)

}

write_rds(error_mnist_tsne, "data/mnist/error_mnist_tsne.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_mnist <- 5:66 #sqrt(NROW(training_data_mnist)/2)

## For tsne
tsne_mnist <- read_rds("data/mnist/mnist_tsne89.rds")
mnist_scaled_obj_tsne <- gen_scaled_data(
  data = tsne_mnist)
tsne_mnist_scaled <- mnist_scaled_obj_tsne$scaled_nldr

lim1 <- mnist_scaled_obj_tsne$lim1
lim2 <- mnist_scaled_obj_tsne$lim2
r2_tsne <- diff(lim2)/diff(lim1)

error_mnist_tsne <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_mnist) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_tsne, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  mnist_model <- fit_highd_model(
    training_data = training_data_mnist,
    emb_df = tsne_mnist_scaled,
    bin1 = xbins,
    r2 = r2_tsne,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = TRUE,
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
    type_NLDR = "tSNE",
    col_start = "PC") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_mnist),
           method = "tSNE2",
           a1 = a1)

  error_mnist_tsne <- bind_rows(error_mnist_tsne, error_df)

}

write_rds(error_mnist_tsne, "data/mnist/error_mnist_tsne2.rds")


###########

## To initialize number of bins along the x-axis
bin1_vec_mnist <- 5:63 #sqrt(NROW(training_data_mnist)/2)

## For phate
phate_mnist <- read_rds("data/mnist/mnist_phate.rds")
mnist_scaled_obj_phate <- gen_scaled_data(
  data = phate_mnist)
phate_mnist_scaled <- mnist_scaled_obj_phate$scaled_nldr

lim1 <- mnist_scaled_obj_phate$lim1
lim2 <- mnist_scaled_obj_phate$lim2
r2_phate <- diff(lim2)/diff(lim1)

error_mnist_phate <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_mnist) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_phate, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  mnist_model <- fit_highd_model(
    training_data = training_data_mnist,
    emb_df = phate_mnist_scaled,
    bin1 = xbins,
    r2 = r2_phate,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = TRUE,
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
    type_NLDR = "PHATE",
    col_start = "PC") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_mnist),
           method = "PHATE",
           a1 = a1)

  error_mnist_phate <- bind_rows(error_mnist_phate, error_df)

}

write_rds(error_mnist_phate, "data/mnist/error_mnist_phate.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_mnist <- 5:63 #sqrt(NROW(training_data_mnist)/2)

## For trimap
trimap_mnist <- read_rds("data/mnist/mnist_trimap.rds")
mnist_scaled_obj_trimap <- gen_scaled_data(
  data = trimap_mnist)
trimap_mnist_scaled <- mnist_scaled_obj_trimap$scaled_nldr

lim1 <- mnist_scaled_obj_trimap$lim1
lim2 <- mnist_scaled_obj_trimap$lim2
r2_trimap <- diff(lim2)/diff(lim1)

error_mnist_trimap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_mnist) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_trimap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1


  mnist_model <- fit_highd_model(
    training_data = training_data_mnist,
    emb_df = trimap_mnist_scaled,
    bin1 = xbins,
    r2 = r2_trimap,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = TRUE,
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
    type_NLDR = "TriMAP",
    col_start = "PC") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_mnist),
           method = "TriMAP",
           a1 = a1)

  error_mnist_trimap <- bind_rows(error_mnist_trimap, error_df)

}

write_rds(error_mnist_trimap, "data/mnist/error_mnist_trimap.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_mnist <- 5:63 #sqrt(NROW(training_data_mnist)/2)

## For pacmap
pacmap_mnist <- read_rds("data/mnist/mnist_pacmap.rds")
mnist_scaled_obj_pacmap <- gen_scaled_data(
  data = pacmap_mnist)
pacmap_mnist_scaled <- mnist_scaled_obj_pacmap$scaled_nldr

lim1 <- mnist_scaled_obj_pacmap$lim1
lim2 <- mnist_scaled_obj_pacmap$lim2
r2_pacmap <- diff(lim2)/diff(lim1)

error_mnist_pacmap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_mnist) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_pacmap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  mnist_model <- fit_highd_model(
    training_data = training_data_mnist,
    emb_df = pacmap_mnist_scaled,
    bin1 = xbins,
    r2 = r2_pacmap,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = TRUE,
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
           b_non_empty = NROW(df_bin_centroids_mnist),
           method = "PaCMAP",
           a1 = a1)

  error_mnist_pacmap <- bind_rows(error_mnist_pacmap, error_df)

}

write_rds(error_mnist_pacmap, "data/mnist/error_mnist_pacmap.rds")
