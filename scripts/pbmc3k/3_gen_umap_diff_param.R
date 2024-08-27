library(readr)
library(uwot)
library(dplyr)
library(ggplot2)
library(quollr)
set.seed(20240110)

n_neighbors <- 9
min_dist <- 0.5

## Select PCs
training_data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50.rds")
training_data_pbmc <- training_data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(training_data_pbmc))

## Obtain UMAP
umap_pbmc <- uwot::umap(training_data_pbmc |> select(-ID),
                  n_neighbors = n_neighbors,
                  n_components =  2,
                  metric = "cosine",
                  min_dist = min_dist,
                  init = "spca")

umap_pbmc <- umap_pbmc |>
  as.data.frame()

names(umap_pbmc)[1:(ncol(umap_pbmc))] <- paste0(rep("UMAP",(ncol(umap_pbmc))), 1:(ncol(umap_pbmc)))

umap_pbmc <- umap_pbmc |>
  mutate(ID = 1:NROW(umap_pbmc))

## Run only once
write_rds(umap_pbmc, file = paste0("data/pbmc3k/pbmc_umap_", n_neighbors, "_min_dist_", min_dist, ".rds"))

## Scaled data
pbmc_scaled_obj <- gen_scaled_data(
  data = umap_pbmc)
umap_pbmc_scaled <- pbmc_scaled_obj$scaled_nldr


## To initialise number of bins along the x-axis
bin1_vec <- 2:40

lim1 <- pbmc_scaled_obj$lim1
lim2 <- pbmc_scaled_obj$lim2
r2_pbmc <- diff(lim2)/diff(lim1)

error_pbmc <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec) {

  bin2 <- calc_bins_y(bin1 = xbins, q = 0.1, r2 = r2_pbmc)$bin2

  pbmc_model <- fit_highd_model(
    training_data = training_data_pbmc,
    emb_df = umap_pbmc_scaled,
    bin1 = xbins,
    q = 0.1,
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


error_pbmc <- error_pbmc |>
  mutate(method = "umap",
         n_neighbors = n_neighbors,
         min_dist = min_dist)


ggplot(error_pbmc, aes(x = b_non_empty,
                       y = log(MSE))) +
  geom_point() +
  geom_line()

write_csv(error_pbmc, "data/pbmc3k/error_df_umap.csv", append = TRUE)

