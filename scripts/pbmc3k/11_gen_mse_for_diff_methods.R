library(readr)
library(quollr)
library(dplyr)

training_data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50.rds")
training_data_pbmc <- training_data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(training_data_pbmc))

names(training_data_pbmc) <- append(paste0("x", 1:9), "ID")

## For umap
umap_pbmc <- read_rds("data/pbmc3k/pbmc_umap_30_min_dist_0.3.rds")
names(umap_pbmc) <- c("emb1", "emb2", "ID")

pbmc_scaled_obj_umap <- gen_scaled_data(
  data = umap_pbmc)
umap_pbmc_scaled <- pbmc_scaled_obj_umap$scaled_nldr

lim1 <- pbmc_scaled_obj_umap$lim1
lim2 <- pbmc_scaled_obj_umap$lim2
r2_umap <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:49 #sqrt(NROW(training_data_pbmc)/r2_umap)

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

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_umap_30_min_dist_0.3.rds")

###########

## For umap
umap_pbmc <- read_rds("data/pbmc3k/pbmc_umap_5_min_dist_0.01.rds")
names(umap_pbmc) <- c("emb1", "emb2", "ID")

pbmc_scaled_obj_umap <- gen_scaled_data(
  data = umap_pbmc)
umap_pbmc_scaled <- pbmc_scaled_obj_umap$scaled_nldr

lim1 <- pbmc_scaled_obj_umap$lim1
lim2 <- pbmc_scaled_obj_umap$lim2
r2_umap <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:57 #sqrt(NROW(training_data_pbmc)/r2_umap)

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
           method = "UMAP_5_min_dist_0.01",
           a1 = a1)

  error_pbmc_umap <- bind_rows(error_pbmc_umap, error_df)

}

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_umap_5_min_dist_0.01.rds")

###########

## For umap

umap_pbmc <- read_rds("data/pbmc3k/pbmc_umap_15_min_dist_0.99.rds")
names(umap_pbmc) <- c("emb1", "emb2", "ID")

pbmc_scaled_obj_umap <- gen_scaled_data(
  data = umap_pbmc)
umap_pbmc_scaled <- pbmc_scaled_obj_umap$scaled_nldr

lim1 <- pbmc_scaled_obj_umap$lim1
lim2 <- pbmc_scaled_obj_umap$lim2
r2_umap <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:45 #sqrt(NROW(training_data_pbmc)/r2_umap)

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
           method = "UMAP_15_min_dist_0.99",
           a1 = a1)

  error_pbmc_umap <- bind_rows(error_pbmc_umap, error_df)

}

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_umap_15_min_dist_0.99.rds")

###########

## For tsne
tsne_pbmc <- read_rds("data/pbmc3k/pbmc_tsne_5.rds")
names(tsne_pbmc) <- c("emb1", "emb2", "ID")

pbmc_scaled_obj_tsne <- gen_scaled_data(
  data = tsne_pbmc)
tsne_pbmc_scaled <- pbmc_scaled_obj_tsne$scaled_nldr

lim1 <- pbmc_scaled_obj_tsne$lim1
lim2 <- pbmc_scaled_obj_tsne$lim2
r2_tsne <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:53 #sqrt(NROW(training_data_pbmc)/r2_tsne)

error_pbmc_umap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

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
           method = "tsne_5",
           a1 = a1)

  error_pbmc_umap <- bind_rows(error_pbmc_umap, error_df)

}

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_tsne_5.rds")

###########

## For tsne
tsne_pbmc <- read_rds("data/pbmc3k/pbmc_tsne_30.rds")
names(tsne_pbmc) <- c("emb1", "emb2", "ID")

pbmc_scaled_obj_tsne <- gen_scaled_data(
  data = tsne_pbmc)
tsne_pbmc_scaled <- pbmc_scaled_obj_tsne$scaled_nldr

lim1 <- pbmc_scaled_obj_tsne$lim1
lim2 <- pbmc_scaled_obj_tsne$lim2
r2_tsne <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:48 #sqrt(NROW(training_data_pbmc)/r2_tsne)

error_pbmc_umap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

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
           method = "tsne_30",
           a1 = a1)

  error_pbmc_umap <- bind_rows(error_pbmc_umap, error_df)

}

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_tsne_30.rds")

###########

## For phate
phate_pbmc <- read_rds("data/pbmc3k/pbmc_phate_5.rds")
names(phate_pbmc) <- c("emb1", "emb2", "ID")

pbmc_scaled_obj_phate <- gen_scaled_data(
  data = phate_pbmc)
phate_pbmc_scaled <- pbmc_scaled_obj_phate$scaled_nldr

lim1 <- pbmc_scaled_obj_phate$lim1
lim2 <- pbmc_scaled_obj_phate$lim2
r2_phate <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:75 #sqrt(NROW(training_data_pbmc)/r2_phate)

error_pbmc_umap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    highd_data = training_data_pbmc,
    nldr_data = phate_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_phate,
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
           method = "phate_5",
           a1 = a1)

  error_pbmc_umap <- bind_rows(error_pbmc_umap, error_df)

}

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_phate_5.rds")

###########

## For trimap
trimap_pbmc <- read_rds("data/pbmc3k/pbmc_trimap_12_4_3.rds")
names(trimap_pbmc) <- c("emb1", "emb2", "ID")

pbmc_scaled_obj_trimap <- gen_scaled_data(
  data = trimap_pbmc)
trimap_pbmc_scaled <- pbmc_scaled_obj_trimap$scaled_nldr

lim1 <- pbmc_scaled_obj_trimap$lim1
lim2 <- pbmc_scaled_obj_trimap$lim2
r2_trimap <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:61 #sqrt(NROW(training_data_pbmc)/r2_trimap)

error_pbmc_umap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    highd_data = training_data_pbmc,
    nldr_data = umap_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_trimap,
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
           method = "trimap_12_4_3",
           a1 = a1)

  error_pbmc_umap <- bind_rows(error_pbmc_umap, error_df)

}

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_trimap_12_4_3.rds")

###########

## For pacmap
pacmap_pbmc <- read_rds("data/pbmc3k/pbmc_pacmap_30_random_0.9_5.rds")
names(pacmap_pbmc) <- c("emb1", "emb2", "ID")

pbmc_scaled_obj_pacmap <- gen_scaled_data(
  data = pacmap_pbmc)
pacmap_pbmc_scaled <- pbmc_scaled_obj_pacmap$scaled_nldr

lim1 <- pbmc_scaled_obj_pacmap$lim1
lim2 <- pbmc_scaled_obj_pacmap$lim2
r2_pacmap <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_pbmc <- 2:58 #sqrt(NROW(training_data_pbmc)/r2_pacmap)

error_pbmc_umap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_pbmc) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  pbmc_model <- fit_highd_model(
    highd_data = training_data_pbmc,
    nldr_data = pacmap_pbmc_scaled,
    bin1 = xbins,
    r2 = r2_pacmap,
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
           method = "pacmap",
           a1 = a1)

  error_pbmc_umap <- bind_rows(error_pbmc_umap, error_df)

}

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_pacmap_30_random_0.9_5.rds")

###########

## For umap
umap_pbmc <- read_rds("data/pbmc3k/pbmc_umap_9_min_dist_0.5.rds")
names(umap_pbmc) <- c("emb1", "emb2", "ID")

pbmc_scaled_obj_umap <- gen_scaled_data(
  data = umap_pbmc)
umap_pbmc_scaled <- pbmc_scaled_obj_umap$scaled_nldr

lim1 <- pbmc_scaled_obj_umap$lim1
lim2 <- pbmc_scaled_obj_umap$lim2
r2_umap <- diff(lim2)/diff(lim1)

error_pbmc_umap <- data.frame(matrix(nrow = 0, ncol = 0))

bin1_vec_pbmc <- 2:50 #sqrt(NROW(training_data_pbmc)/r2_umap)

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
           method = "UMAP_9_min_dist_0.5",
           a1 = a1)

  error_pbmc_umap <- bind_rows(error_pbmc_umap, error_df)

}

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_umap_9_min_dist_0.5.rds")
