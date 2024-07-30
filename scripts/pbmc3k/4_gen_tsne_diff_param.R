library(readr)
library(Rtsne)
library(dplyr)
library(ggplot2)
library(quollr)
set.seed(20240110)
source("nldr_code.R")

perplexity <- 51

## Select PCs
training_data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50.rds")
training_data_pbmc <- training_data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(training_data_pbmc))

## Obtain tSNE
tSNE_data <- Fit_tSNE(training_data_pbmc |> select(-ID),
                      opt_perplexity = perplexity,
                      with_seed = 20240110)

tSNE_data <- tSNE_data |>
  select(-ID) |>
  mutate(ID = training_data_pbmc$ID)

## Run only once
write_rds(tSNE_data, file = paste0("data/pbmc3k/pbmc_tsne_", perplexity, ".rds"))


## Scaled data
pbmc_scaled_obj <- gen_scaled_data(
  data = tSNE_data)
tsne_pbmc_scaled <- pbmc_scaled_obj$scaled_nldr


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
    emb_df = tsne_pbmc_scaled,
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
    type_NLDR = "tSNE",
    col_start = "PC_") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_pbmc))

  error_pbmc <- bind_rows(error_pbmc, error_df)

}


error_pbmc <- error_pbmc |>
  mutate(method = "tsne",
         perplexity = perplexity)


ggplot(error_pbmc, aes(x = b_non_empty,
                       y = log(MSE))) +
  geom_point() +
  geom_line()

write_csv(error_pbmc, "data/pbmc3k/error_df_tsne.csv", append = TRUE)
