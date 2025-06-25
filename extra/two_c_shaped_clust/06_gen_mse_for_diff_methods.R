library(readr)
library(quollr)
library(dplyr)

quad <- function(a = 3, b = 2 * a2, c = -(a2^2 + a1^2))
{
  a <- as.complex(a)
  answer <- c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
              (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
  if(all(Im(answer) == 0)) answer <- Re(answer)
  if(answer[1] == answer[2]) return(answer[1])
  answer[answer>0] ## only positive
}


training_data_two_curvy_diff_clust <- read_rds("data/two_curvy_diff_clust/two_curvy_diff_clust_data.rds")
training_data_two_curvy_diff_clust <- training_data_two_curvy_diff_clust |>
  mutate(ID = 1:NROW(training_data_two_curvy_diff_clust))

## To initialize number of bins along the x-axis
bin1_vec_two_curvy_diff_clust <- 4:67 #sqrt(NROW(training_data_two_curvy_diff_clust)/r2_umap)

## For umap
umap_two_curvy_diff_clust <- read_rds("data/two_curvy_diff_clust/two_curvy_diff_clust_umap_n-neigbors_15_min-dist_0.1.rds")
two_curvy_diff_clust_scaled_obj_umap <- gen_scaled_data(
  data = umap_two_curvy_diff_clust)
umap_two_curvy_diff_clust_scaled <- two_curvy_diff_clust_scaled_obj_umap$scaled_nldr |>
  mutate(ID = 1:NROW(training_data_two_curvy_diff_clust))

lim1 <- two_curvy_diff_clust_scaled_obj_umap$lim1
lim2 <- two_curvy_diff_clust_scaled_obj_umap$lim2
r2_umap <- diff(lim2)/diff(lim1)

error_two_curvy_diff_clust_umap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_two_curvy_diff_clust) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1
  a2 <- hb_obj$a2

  two_curvy_diff_clust_model <- fit_highd_model(
    training_data = training_data_two_curvy_diff_clust,
    emb_df = umap_two_curvy_diff_clust_scaled,
    bin1 = xbins,
    r2 = r2_umap,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "x"
  )

  df_bin_centroids_two_curvy_diff_clust <- two_curvy_diff_clust_model$df_bin_centroids
  df_bin_two_curvy_diff_clust <- two_curvy_diff_clust_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_two_curvy_diff_clust,
    df_bin = df_bin_two_curvy_diff_clust,
    training_data = training_data_two_curvy_diff_clust,
    newdata = NULL,
    type_NLDR = "UMAP",
    col_start = "x") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_two_curvy_diff_clust),
           method = "UMAP",
           a1 = round(a1, 2),
           a2 = round(a2, 2),
           side_length = quad(a=3, b = 2 * a2, c = -(a2^2 + a1^2)))

  error_two_curvy_diff_clust_umap <- bind_rows(error_two_curvy_diff_clust_umap, error_df)

}

write_rds(error_two_curvy_diff_clust_umap, "data/two_curvy_diff_clust/error_two_curvy_diff_clust_umap.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_two_curvy_diff_clust <- 5:67 #sqrt(NROW(training_data_two_curvy_diff_clust)/r2_tsne)

## For tsne
#tsne_two_curvy_diff_clust <- read_rds("data/two_curvy_diff_clust/two_curvy_diff_clust_tsne89.rds")
tsne_two_curvy_diff_clust <- read_rds("data/two_curvy_diff_clust/two_curvy_diff_clust_tsne_perplexity_30.rds")
two_curvy_diff_clust_scaled_obj_tsne <- gen_scaled_data(
  data = tsne_two_curvy_diff_clust)
tsne_two_curvy_diff_clust_scaled <- two_curvy_diff_clust_scaled_obj_tsne$scaled_nldr |>
  mutate(ID = 1:NROW(training_data_two_curvy_diff_clust))

lim1 <- two_curvy_diff_clust_scaled_obj_tsne$lim1
lim2 <- two_curvy_diff_clust_scaled_obj_tsne$lim2
r2_tsne <- diff(lim2)/diff(lim1)

error_two_curvy_diff_clust_tsne <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_two_curvy_diff_clust) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_tsne, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1
  a2 <- hb_obj$a2

  two_curvy_diff_clust_model <- fit_highd_model(
    training_data = training_data_two_curvy_diff_clust,
    emb_df = tsne_two_curvy_diff_clust_scaled,
    bin1 = xbins,
    r2 = r2_tsne,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "x"
  )

  df_bin_centroids_two_curvy_diff_clust <- two_curvy_diff_clust_model$df_bin_centroids
  df_bin_two_curvy_diff_clust <- two_curvy_diff_clust_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_two_curvy_diff_clust,
    df_bin = df_bin_two_curvy_diff_clust,
    training_data = training_data_two_curvy_diff_clust,
    newdata = NULL,
    type_NLDR = "tSNE",
    col_start = "x") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_two_curvy_diff_clust),
           method = "tSNE",
           a1 = round(a1, 2),
           a2 = round(a2, 2),
           side_length = quad(a=3, b = 2 * a2, c = -(a2^2 + a1^2)))

  error_two_curvy_diff_clust_tsne <- bind_rows(error_two_curvy_diff_clust_tsne, error_df)

}

write_rds(error_two_curvy_diff_clust_tsne, "data/two_curvy_diff_clust/error_two_curvy_diff_clust_tsne.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_two_curvy_diff_clust <- 5:67 #sqrt(NROW(training_data_two_curvy_diff_clust)/2)

## For phate
phate_two_curvy_diff_clust <- read_rds("data/two_curvy_diff_clust/two_curvy_diff_clust_phate_knn_5.rds")
two_curvy_diff_clust_scaled_obj_phate <- gen_scaled_data(
  data = phate_two_curvy_diff_clust)
phate_two_curvy_diff_clust_scaled <- two_curvy_diff_clust_scaled_obj_phate$scaled_nldr |>
  mutate(ID = 1:NROW(training_data_two_curvy_diff_clust))

lim1 <- two_curvy_diff_clust_scaled_obj_phate$lim1
lim2 <- two_curvy_diff_clust_scaled_obj_phate$lim2
r2_phate <- diff(lim2)/diff(lim1)

error_two_curvy_diff_clust_phate <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_two_curvy_diff_clust) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_phate, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1
  a2 <- hb_obj$a2

  two_curvy_diff_clust_model <- fit_highd_model(
    training_data = training_data_two_curvy_diff_clust,
    emb_df = phate_two_curvy_diff_clust_scaled,
    bin1 = xbins,
    r2 = r2_phate,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "x"
  )

  df_bin_centroids_two_curvy_diff_clust <- two_curvy_diff_clust_model$df_bin_centroids
  df_bin_two_curvy_diff_clust <- two_curvy_diff_clust_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_two_curvy_diff_clust,
    df_bin = df_bin_two_curvy_diff_clust,
    training_data = training_data_two_curvy_diff_clust,
    newdata = NULL,
    type_NLDR = "PHATE",
    col_start = "x") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_two_curvy_diff_clust),
           method = "PHATE",
           a1 = round(a1, 2),
           a2 = round(a2, 2),
           side_length = quad(a=3, b = 2 * a2, c = -(a2^2 + a1^2)))

  error_two_curvy_diff_clust_phate <- bind_rows(error_two_curvy_diff_clust_phate, error_df)

}

write_rds(error_two_curvy_diff_clust_phate, "data/two_curvy_diff_clust/error_two_curvy_diff_clust_phate.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_two_curvy_diff_clust <- 5:67 #sqrt(NROW(training_data_two_curvy_diff_clust)/2)

## For trimap
trimap_two_curvy_diff_clust <- read_rds("data/two_curvy_diff_clust/two_curvy_diff_clust_trimap_n-inliers_12_n-outliers_4_n-random_3.rds")
two_curvy_diff_clust_scaled_obj_trimap <- gen_scaled_data(
  data = trimap_two_curvy_diff_clust)
trimap_two_curvy_diff_clust_scaled <- two_curvy_diff_clust_scaled_obj_trimap$scaled_nldr |>
  mutate(ID = 1:NROW(training_data_two_curvy_diff_clust))

lim1 <- two_curvy_diff_clust_scaled_obj_trimap$lim1
lim2 <- two_curvy_diff_clust_scaled_obj_trimap$lim2
r2_trimap <- diff(lim2)/diff(lim1)

error_two_curvy_diff_clust_trimap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_two_curvy_diff_clust) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_trimap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1
  a2 <- hb_obj$a2

  two_curvy_diff_clust_model <- fit_highd_model(
    training_data = training_data_two_curvy_diff_clust,
    emb_df = trimap_two_curvy_diff_clust_scaled,
    bin1 = xbins,
    r2 = r2_trimap,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "x"
  )

  df_bin_centroids_two_curvy_diff_clust <- two_curvy_diff_clust_model$df_bin_centroids
  df_bin_two_curvy_diff_clust <- two_curvy_diff_clust_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_two_curvy_diff_clust,
    df_bin = df_bin_two_curvy_diff_clust,
    training_data = training_data_two_curvy_diff_clust,
    newdata = NULL,
    type_NLDR = "TriMAP",
    col_start = "x") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_two_curvy_diff_clust),
           method = "TriMAP",
           a1 = round(a1, 2),
           a2 = round(a2, 2),
           side_length = quad(a=3, b = 2 * a2, c = -(a2^2 + a1^2)))

  error_two_curvy_diff_clust_trimap <- bind_rows(error_two_curvy_diff_clust_trimap, error_df)

}

write_rds(error_two_curvy_diff_clust_trimap, "data/two_curvy_diff_clust/error_two_curvy_diff_clust_trimap.rds")

###########

## To initialize number of bins along the x-axis
bin1_vec_two_curvy_diff_clust <- 5:67 #sqrt(NROW(training_data_two_curvy_diff_clust)/2)

## For pacmap
pacmap_two_curvy_diff_clust <- read_rds("data/two_curvy_diff_clust/two_curvy_diff_clust_pacmap_n-neighbors_10_init_random_MN-ratio_0.5_FP-ratio_2.rds")
two_curvy_diff_clust_scaled_obj_pacmap <- gen_scaled_data(
  data = pacmap_two_curvy_diff_clust)
pacmap_two_curvy_diff_clust_scaled <- two_curvy_diff_clust_scaled_obj_pacmap$scaled_nldr |>
  mutate(ID = 1:NROW(training_data_two_curvy_diff_clust))

lim1 <- two_curvy_diff_clust_scaled_obj_pacmap$lim1
lim2 <- two_curvy_diff_clust_scaled_obj_pacmap$lim2
r2_pacmap <- diff(lim2)/diff(lim1)

error_two_curvy_diff_clust_pacmap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_two_curvy_diff_clust) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_pacmap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1
  a2 <- hb_obj$a2

  two_curvy_diff_clust_model <- fit_highd_model(
    training_data = training_data_two_curvy_diff_clust,
    emb_df = pacmap_two_curvy_diff_clust_scaled,
    bin1 = xbins,
    r2 = r2_pacmap,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "x"
  )

  df_bin_centroids_two_curvy_diff_clust <- two_curvy_diff_clust_model$df_bin_centroids
  df_bin_two_curvy_diff_clust <- two_curvy_diff_clust_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_two_curvy_diff_clust,
    df_bin = df_bin_two_curvy_diff_clust,
    training_data = training_data_two_curvy_diff_clust,
    newdata = NULL,
    type_NLDR = "PaCMAP",
    col_start = "x") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_two_curvy_diff_clust),
           method = "PaCMAP",
           a1 = round(a1, 2),
           a2 = round(a2, 2),
           side_length = quad(a=3, b = 2 * a2, c = -(a2^2 + a1^2)))

  error_two_curvy_diff_clust_pacmap <- bind_rows(error_two_curvy_diff_clust_pacmap, error_df)

}

write_rds(error_two_curvy_diff_clust_pacmap, "data/two_curvy_diff_clust/error_two_curvy_diff_clust_pacmap.rds")
