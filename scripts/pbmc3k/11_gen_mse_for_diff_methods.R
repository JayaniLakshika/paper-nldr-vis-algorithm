library(readr)
library(quollr)
library(dplyr)

training_data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50.rds")
training_data_pbmc <- training_data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(training_data_pbmc))

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:36 #sqrt(NROW(training_data_pbmc)/2)

## For umap
umap_pbmc <- read_rds("data/pbmc3k/pbmc_umap_30_min_dist_0.3.rds")
pbmc_scaled_obj_umap <- gen_scaled_data(
  data = umap_pbmc)
umap_pbmc_scaled <- pbmc_scaled_obj_umap$scaled_nldr

lim1 <- pbmc_scaled_obj_umap$lim1
lim2 <- pbmc_scaled_obj_umap$lim2
r2_umap <- diff(lim2)/diff(lim1)

error_pbmc_umap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    training_data = training_data_pbmc,
    emb_df = umap_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_umap,
    q = 0.1,
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
           b_non_empty = NROW(df_bin_centroids_pbmc),
           method = "UMAP_30_min_dist_0.3",
           a1 = a1)

  error_pbmc_umap <- bind_rows(error_pbmc_umap, error_df)

}

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_umap_30_min_dist_0.3.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:36 #sqrt(NROW(training_data_pbmc)/2)

## For umap
umap_pbmc <- read_rds("data/pbmc3k/pbmc_umap_5_min_dist_0.01.rds")
pbmc_scaled_obj_umap <- gen_scaled_data(
  data = umap_pbmc)
umap_pbmc_scaled <- pbmc_scaled_obj_umap$scaled_nldr

lim1 <- pbmc_scaled_obj_umap$lim1
lim2 <- pbmc_scaled_obj_umap$lim2
r2_umap <- diff(lim2)/diff(lim1)

error_pbmc_umap2 <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    training_data = training_data_pbmc,
    emb_df = umap_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_umap,
    q = 0.1,
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
           b_non_empty = NROW(df_bin_centroids_pbmc),
           method = "UMAP_5_min_dist_0.01",
           a1 = a1)

  error_pbmc_umap2 <- bind_rows(error_pbmc_umap2, error_df)

}

write_rds(error_pbmc_umap2, "data/pbmc3k/error_pbmc_umap_5_min_dist_0.01.rds")

###########

## For umap

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:36 #sqrt(NROW(training_data_pbmc)/2)
umap_pbmc <- read_rds("data/pbmc3k/pbmc_umap_15_min_dist_0.99.rds")
pbmc_scaled_obj_umap <- gen_scaled_data(
  data = umap_pbmc)
umap_pbmc_scaled <- pbmc_scaled_obj_umap$scaled_nldr

lim1 <- pbmc_scaled_obj_umap$lim1
lim2 <- pbmc_scaled_obj_umap$lim2
r2_umap <- diff(lim2)/diff(lim1)

error_pbmc_umap3 <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    training_data = training_data_pbmc,
    emb_df = umap_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_umap,
    q = 0.1,
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
           b_non_empty = NROW(df_bin_centroids_pbmc),
           method = "UMAP_15_min_dist_0.99",
           a1 = a1)

  error_pbmc_umap3 <- bind_rows(error_pbmc_umap3, error_df)

}

write_rds(error_pbmc_umap3, "data/pbmc3k/error_pbmc_umap_15_min_dist_0.99.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:22 #sqrt(NROW(training_data_pbmc)/2)

## For tsne
tsne_pbmc <- read_rds("data/pbmc3k/pbmc_tsne_5.rds")
pbmc_scaled_obj_tsne <- gen_scaled_data(
  data = tsne_pbmc)
tsne_pbmc_scaled <- pbmc_scaled_obj_tsne$scaled_nldr

lim1 <- pbmc_scaled_obj_tsne$lim1
lim2 <- pbmc_scaled_obj_tsne$lim2
r2_tsne <- diff(lim2)/diff(lim1)

error_pbmc_tsne <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_tsne, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    training_data = training_data_pbmc,
    emb_df = tsne_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_tsne,
    q = 0.1,
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
    type_NLDR = "tSNE",
    col_start = "PC_") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_pbmc),
           method = "tsne_5",
           a1 = a1)

  error_pbmc_tsne <- bind_rows(error_pbmc_tsne, error_df)

}

write_rds(error_pbmc_tsne, "data/pbmc3k/error_pbmc_tsne_5.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:36 #sqrt(NROW(training_data_pbmc)/2)

## For tsne
tsne_pbmc <- read_rds("data/pbmc3k/pbmc_tsne_30.rds")
pbmc_scaled_obj_tsne <- gen_scaled_data(
  data = tsne_pbmc)
tsne_pbmc_scaled <- pbmc_scaled_obj_tsne$scaled_nldr

lim1 <- pbmc_scaled_obj_tsne$lim1
lim2 <- pbmc_scaled_obj_tsne$lim2
r2_tsne <- diff(lim2)/diff(lim1)

error_pbmc_tsne2 <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_tsne, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    training_data = training_data_pbmc,
    emb_df = tsne_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_tsne,
    q = 0.1,
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
    type_NLDR = "tSNE",
    col_start = "PC_") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_pbmc),
           method = "tsne_30",
           a1 = a1)

  error_pbmc_tsne2 <- bind_rows(error_pbmc_tsne2, error_df)

}

write_rds(error_pbmc_tsne2, "data/pbmc3k/error_pbmc_tsne_30.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:46 #sqrt(NROW(training_data_pbmc)/2)

## For phate
phate_pbmc <- read_rds("data/pbmc3k/pbmc_phate_5.rds")
pbmc_scaled_obj_phate <- gen_scaled_data(
  data = phate_pbmc)
phate_pbmc_scaled <- pbmc_scaled_obj_phate$scaled_nldr

lim1 <- pbmc_scaled_obj_phate$lim1
lim2 <- pbmc_scaled_obj_phate$lim2
r2_phate <- diff(lim2)/diff(lim1)

error_pbmc_phate <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_phate, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    training_data = training_data_pbmc,
    emb_df = phate_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_phate,
    q = 0.1,
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
    type_NLDR = "PHATE",
    col_start = "PC_") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_pbmc),
           method = "phate_5",
           a1 = a1)

  error_pbmc_phate <- bind_rows(error_pbmc_phate, error_df)

}

write_rds(error_pbmc_phate, "data/pbmc3k/error_pbmc_phate_5.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:36 #sqrt(NROW(training_data_pbmc)/2)

## For trimap
trimap_pbmc <- read_rds("data/pbmc3k/pbmc_trimap_12_4_3.rds")
pbmc_scaled_obj_trimap <- gen_scaled_data(
  data = trimap_pbmc)
trimap_pbmc_scaled <- pbmc_scaled_obj_trimap$scaled_nldr

lim1 <- pbmc_scaled_obj_trimap$lim1
lim2 <- pbmc_scaled_obj_trimap$lim2
r2_trimap <- diff(lim2)/diff(lim1)

error_pbmc_trimap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_trimap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1


  pbmc_model <- fit_highd_model(
    training_data = training_data_pbmc,
    emb_df = trimap_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_trimap,
    q = 0.1,
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
    type_NLDR = "TriMAP",
    col_start = "PC_") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_pbmc),
           method = "trimap_12_4_3",
           a1 = a1)

  error_pbmc_trimap <- bind_rows(error_pbmc_trimap, error_df)

}

write_rds(error_pbmc_trimap, "data/pbmc3k/error_pbmc_trimap_12_4_3.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:43 #sqrt(NROW(training_data_pbmc)/2)

## For pacmap
pacmap_pbmc <- read_rds("data/pbmc3k/pbmc_pacmap_30_random_0.9_5.rds")
pbmc_scaled_obj_pacmap <- gen_scaled_data(
  data = pacmap_pbmc)
pacmap_pbmc_scaled <- pbmc_scaled_obj_pacmap$scaled_nldr

lim1 <- pbmc_scaled_obj_pacmap$lim1
lim2 <- pbmc_scaled_obj_pacmap$lim2
r2_pacmap <- diff(lim2)/diff(lim1)

error_pbmc_pacmap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_pacmap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    training_data = training_data_pbmc,
    emb_df = pacmap_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_pacmap,
    q = 0.1,
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
    type_NLDR = "PaCMAP",
    col_start = "PC_") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_pbmc),
           method = "pacmap",
           a1 = a1)

  error_pbmc_pacmap <- bind_rows(error_pbmc_pacmap, error_df)

}

write_rds(error_pbmc_pacmap, "data/pbmc3k/error_pbmc_pacmap_30_random_0.9_5.rds")

###########
