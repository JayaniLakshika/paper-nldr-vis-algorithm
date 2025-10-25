## This script is to generate MSE for hyper-parameter in UMAP suggested by Chen and scDEED
library(readr)
library(quollr)
library(dplyr)

data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50_scdeed.rds")
data_pbmc <- data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(data_pbmc))

names(data_pbmc) <- append(paste0("x", 1:9), "ID")

## For umap
umap_pbmc <- read_rds("data/pbmc3k/pbmc_scdeed_umap_n_neighbors_30_min_dist_0.3.rds")
umap_pbmc <- as_tibble(umap_pbmc)
names(umap_pbmc) <- c("emb1", "emb2")
umap_pbmc <- umap_pbmc |>
  mutate(ID = 1:NROW(umap_pbmc))

error_pbmc_umap <- gen_diffbin1_errors(highd_data = data_pbmc, nldr_data = umap_pbmc) |>
  dplyr::mutate(method = "UMAP_30_min_dist_0.3")

write_rds(error_pbmc_umap, "data/pbmc3k/error_scdeed_pbmc_umap_30_min_dist_0.3.rds")

#################################

## For umap
umap_pbmc <- read_rds("data/pbmc3k/pbmc_scdeed_umap_n_neighbors_80_min_dist_0.5.rds")
umap_pbmc <- as_tibble(umap_pbmc)
names(umap_pbmc) <- c("emb1", "emb2")
umap_pbmc <- umap_pbmc |>
  mutate(ID = 1:NROW(umap_pbmc))

error_pbmc_umap <- gen_diffbin1_errors(highd_data = data_pbmc, nldr_data = umap_pbmc) |>
  dplyr::mutate(method = "UMAP_80_min_dist_0.5")

write_rds(error_pbmc_umap, "data/pbmc3k/error_scdeed_pbmc_umap_80_min_dist_0.5.rds")

