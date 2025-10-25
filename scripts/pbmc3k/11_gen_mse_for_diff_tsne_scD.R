## This script is to generate MSE for hyper-parameters in tSNE suggested by Chen and scDEED

library(readr)
library(quollr)
library(dplyr)

data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50_scdeed.rds")
data_pbmc <- data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(data_pbmc))

names(data_pbmc) <- append(paste0("x", 1:9), "ID")

## For tsne
tsne_pbmc <- read_rds("data/pbmc3k/pbmc_scdeed_tsne_perplexity_30.rds")
tsne_pbmc <- as_tibble(tsne_pbmc)
names(tsne_pbmc) <- c("emb1", "emb2")
tsne_pbmc <- tsne_pbmc |>
  mutate(ID = 1:NROW(tsne_pbmc))

error_pbmc_tsne <- gen_diffbin1_errors(highd_data = data_pbmc, nldr_data = tsne_pbmc) |>
  dplyr::mutate(method = "tsne_perplexity_30")

write_rds(error_pbmc_tsne, "data/pbmc3k/error_scdeed_pbmc_tsne_perplexity_30.rds")

#################################

## For tsne
tsne_pbmc <- read_rds("data/pbmc3k/pbmc_scdeed_tsne_perplexity_320.rds")
tsne_pbmc <- as_tibble(tsne_pbmc)
names(tsne_pbmc) <- c("emb1", "emb2")
tsne_pbmc <- tsne_pbmc |>
  mutate(ID = 1:NROW(tsne_pbmc))

error_pbmc_tsne <- gen_diffbin1_errors(highd_data = data_pbmc, nldr_data = tsne_pbmc) |>
  dplyr::mutate(method = "tsne_perplexity_320")

write_rds(error_pbmc_tsne, "data/pbmc3k/error_scdeed_pbmc_tsne_perplexity_320.rds")

