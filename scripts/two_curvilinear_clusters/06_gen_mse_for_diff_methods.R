library(readr)
library(quollr)
library(dplyr)

training_data_two_nonlinear_clusters <- read_rds("data/two_nonlinear_clusters/two_nonlinear_clusters_data.rds")
training_data_two_nonlinear_clusters <- training_data_two_nonlinear_clusters |>
  mutate(ID = 1:NROW(training_data_two_nonlinear_clusters))

## To initialize number of bins along the x-axis
bin1_vec_two_nonlinear_clusters <- 5:22 #sqrt(NROW(training_data_two_nonlinear_clusters)/2)

## For umap
umap_two_nonlinear_clusters <- read_rds("data/two_nonlinear_clusters/two_nonlinear_clusters_umap_n-neigbors_15_min-dist_0.1.rds")
two_nonlinear_clusters_scaled_obj_umap <- gen_scaled_data(
  data = umap_two_nonlinear_clusters)
umap_two_nonlinear_clusters_scaled <- two_nonlinear_clusters_scaled_obj_umap$scaled_nldr |>
  mutate(ID = 1:NROW(training_data_two_nonlinear_clusters))

lim1 <- two_nonlinear_clusters_scaled_obj_umap$lim1
lim2 <- two_nonlinear_clusters_scaled_obj_umap$lim2
r2_umap <- diff(lim2)/diff(lim1)

error_two_nonlinear_clusters_umap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_two_nonlinear_clusters) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  two_nonlinear_clusters_model <- fit_highd_model(
    training_data = training_data_two_nonlinear_clusters,
    emb_df = umap_two_nonlinear_clusters_scaled,
    bin1 = xbins,
    r2 = r2_umap,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "x"
  )

  df_bin_centroids_two_nonlinear_clusters <- two_nonlinear_clusters_model$df_bin_centroids
  df_bin_two_nonlinear_clusters <- two_nonlinear_clusters_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_two_nonlinear_clusters,
    df_bin = df_bin_two_nonlinear_clusters,
    training_data = training_data_two_nonlinear_clusters,
    newdata = NULL,
    type_NLDR = "UMAP",
    col_start = "x") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_two_nonlinear_clusters),
           method = "UMAP",
           a1 = a1)

  error_two_nonlinear_clusters_umap <- bind_rows(error_two_nonlinear_clusters_umap, error_df)

}

write_rds(error_two_nonlinear_clusters_umap, "data/two_nonlinear_clusters/error_two_nonlinear_clusters_umap.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_two_nonlinear_clusters <- 5:22 #sqrt(NROW(training_data_two_nonlinear_clusters)/2)

## For tsne
#tsne_two_nonlinear_clusters <- read_rds("data/two_nonlinear_clusters/two_nonlinear_clusters_tsne89.rds")
tsne_two_nonlinear_clusters <- read_rds("data/two_nonlinear_clusters/two_nonlinear_clusters_tsne_perplexity_30.rds")
two_nonlinear_clusters_scaled_obj_tsne <- gen_scaled_data(
  data = tsne_two_nonlinear_clusters)
tsne_two_nonlinear_clusters_scaled <- two_nonlinear_clusters_scaled_obj_tsne$scaled_nldr |>
  mutate(ID = 1:NROW(training_data_two_nonlinear_clusters))

lim1 <- two_nonlinear_clusters_scaled_obj_tsne$lim1
lim2 <- two_nonlinear_clusters_scaled_obj_tsne$lim2
r2_tsne <- diff(lim2)/diff(lim1)

error_two_nonlinear_clusters_tsne <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_two_nonlinear_clusters) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_tsne, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  two_nonlinear_clusters_model <- fit_highd_model(
    training_data = training_data_two_nonlinear_clusters,
    emb_df = tsne_two_nonlinear_clusters_scaled,
    bin1 = xbins,
    r2 = r2_tsne,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "x"
  )

  df_bin_centroids_two_nonlinear_clusters <- two_nonlinear_clusters_model$df_bin_centroids
  df_bin_two_nonlinear_clusters <- two_nonlinear_clusters_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_two_nonlinear_clusters,
    df_bin = df_bin_two_nonlinear_clusters,
    training_data = training_data_two_nonlinear_clusters,
    newdata = NULL,
    type_NLDR = "tSNE",
    col_start = "x") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_two_nonlinear_clusters),
           method = "tSNE",
           a1 = a1)

  error_two_nonlinear_clusters_tsne <- bind_rows(error_two_nonlinear_clusters_tsne, error_df)

}

write_rds(error_two_nonlinear_clusters_tsne, "data/two_nonlinear_clusters/error_two_nonlinear_clusters_tsne.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_two_nonlinear_clusters <- 5:22 #sqrt(NROW(training_data_two_nonlinear_clusters)/2)

## For phate
phate_two_nonlinear_clusters <- read_rds("data/two_nonlinear_clusters/two_nonlinear_clusters_phate_knn_5.rds")
two_nonlinear_clusters_scaled_obj_phate <- gen_scaled_data(
  data = phate_two_nonlinear_clusters)
phate_two_nonlinear_clusters_scaled <- two_nonlinear_clusters_scaled_obj_phate$scaled_nldr |>
  mutate(ID = 1:NROW(training_data_two_nonlinear_clusters))

lim1 <- two_nonlinear_clusters_scaled_obj_phate$lim1
lim2 <- two_nonlinear_clusters_scaled_obj_phate$lim2
r2_phate <- diff(lim2)/diff(lim1)

error_two_nonlinear_clusters_phate <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_two_nonlinear_clusters) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_phate, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  two_nonlinear_clusters_model <- fit_highd_model(
    training_data = training_data_two_nonlinear_clusters,
    emb_df = phate_two_nonlinear_clusters_scaled,
    bin1 = xbins,
    r2 = r2_phate,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "x"
  )

  df_bin_centroids_two_nonlinear_clusters <- two_nonlinear_clusters_model$df_bin_centroids
  df_bin_two_nonlinear_clusters <- two_nonlinear_clusters_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_two_nonlinear_clusters,
    df_bin = df_bin_two_nonlinear_clusters,
    training_data = training_data_two_nonlinear_clusters,
    newdata = NULL,
    type_NLDR = "PHATE",
    col_start = "x") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_two_nonlinear_clusters),
           method = "PHATE",
           a1 = a1)

  error_two_nonlinear_clusters_phate <- bind_rows(error_two_nonlinear_clusters_phate, error_df)

}

write_rds(error_two_nonlinear_clusters_phate, "data/two_nonlinear_clusters/error_two_nonlinear_clusters_phate.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_two_nonlinear_clusters <- 5:22 #sqrt(NROW(training_data_two_nonlinear_clusters)/2)

## For trimap
trimap_two_nonlinear_clusters <- read_rds("data/two_nonlinear_clusters/two_nonlinear_clusters_trimap_n-inliers_12_n-outliers_4_n-random_3.rds")
two_nonlinear_clusters_scaled_obj_trimap <- gen_scaled_data(
  data = trimap_two_nonlinear_clusters)
trimap_two_nonlinear_clusters_scaled <- two_nonlinear_clusters_scaled_obj_trimap$scaled_nldr |>
  mutate(ID = 1:NROW(training_data_two_nonlinear_clusters))

lim1 <- two_nonlinear_clusters_scaled_obj_trimap$lim1
lim2 <- two_nonlinear_clusters_scaled_obj_trimap$lim2
r2_trimap <- diff(lim2)/diff(lim1)

error_two_nonlinear_clusters_trimap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_two_nonlinear_clusters) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_trimap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1


  two_nonlinear_clusters_model <- fit_highd_model(
    training_data = training_data_two_nonlinear_clusters,
    emb_df = trimap_two_nonlinear_clusters_scaled,
    bin1 = xbins,
    r2 = r2_trimap,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "x"
  )

  df_bin_centroids_two_nonlinear_clusters <- two_nonlinear_clusters_model$df_bin_centroids
  df_bin_two_nonlinear_clusters <- two_nonlinear_clusters_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_two_nonlinear_clusters,
    df_bin = df_bin_two_nonlinear_clusters,
    training_data = training_data_two_nonlinear_clusters,
    newdata = NULL,
    type_NLDR = "TriMAP",
    col_start = "x") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_two_nonlinear_clusters),
           method = "TriMAP",
           a1 = a1)

  error_two_nonlinear_clusters_trimap <- bind_rows(error_two_nonlinear_clusters_trimap, error_df)

}

write_rds(error_two_nonlinear_clusters_trimap, "data/two_nonlinear_clusters/error_two_nonlinear_clusters_trimap.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_two_nonlinear_clusters <- 5:22 #sqrt(NROW(training_data_two_nonlinear_clusters)/2)

## For pacmap
pacmap_two_nonlinear_clusters <- read_rds("data/two_nonlinear_clusters/two_nonlinear_clusters_pacmap_n-neighbors_10_init_random_MN-ratio_0.5_FP-ratio_2.rds")
two_nonlinear_clusters_scaled_obj_pacmap <- gen_scaled_data(
  data = pacmap_two_nonlinear_clusters)
pacmap_two_nonlinear_clusters_scaled <- two_nonlinear_clusters_scaled_obj_pacmap$scaled_nldr |>
  mutate(ID = 1:NROW(training_data_two_nonlinear_clusters))

lim1 <- two_nonlinear_clusters_scaled_obj_pacmap$lim1
lim2 <- two_nonlinear_clusters_scaled_obj_pacmap$lim2
r2_pacmap <- diff(lim2)/diff(lim1)

error_two_nonlinear_clusters_pacmap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_two_nonlinear_clusters) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_pacmap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  two_nonlinear_clusters_model <- fit_highd_model(
    training_data = training_data_two_nonlinear_clusters,
    emb_df = pacmap_two_nonlinear_clusters_scaled,
    bin1 = xbins,
    r2 = r2_pacmap,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "x"
  )

  df_bin_centroids_two_nonlinear_clusters <- two_nonlinear_clusters_model$df_bin_centroids
  df_bin_two_nonlinear_clusters <- two_nonlinear_clusters_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_two_nonlinear_clusters,
    df_bin = df_bin_two_nonlinear_clusters,
    training_data = training_data_two_nonlinear_clusters,
    newdata = NULL,
    type_NLDR = "PaCMAP",
    col_start = "x") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_two_nonlinear_clusters),
           method = "PaCMAP",
           a1 = a1)

  error_two_nonlinear_clusters_pacmap <- bind_rows(error_two_nonlinear_clusters_pacmap, error_df)

}

write_rds(error_two_nonlinear_clusters_pacmap, "data/two_nonlinear_clusters/error_two_nonlinear_clusters_pacmap.rds")
