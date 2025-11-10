library(readr)
library(quollr)
library(dplyr)

data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50.rds")
data_pbmc <- data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(data_pbmc))

names(data_pbmc) <- append(paste0("x", 1:9), "ID")

## For umap
umap_pbmc <- read_rds("data/pbmc3k/pbmc_umap_30_min_dist_0.3.rds")
names(umap_pbmc) <- c("emb1", "emb2", "ID")

error_pbmc_umap <- gen_diffbin1_errors(highd_data = data_pbmc, nldr_data = umap_pbmc) |>
  dplyr::mutate(method = "UMAP_30_min_dist_0.3")

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_umap_30_min_dist_0.3.rds")

###########

## For umap
umap_pbmc <- read_rds("data/pbmc3k/pbmc_umap_5_min_dist_0.01.rds")
names(umap_pbmc) <- c("emb1", "emb2", "ID")

error_pbmc_umap <- gen_diffbin1_errors(highd_data = data_pbmc, nldr_data = umap_pbmc) |>
  dplyr::mutate(method = "UMAP_5_min_dist_0.01")

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_umap_5_min_dist_0.01.rds")

###########

## For umap

umap_pbmc <- read_rds("data/pbmc3k/pbmc_umap_5_min_dist_0.8.rds")
names(umap_pbmc) <- c("emb1", "emb2", "ID")

error_pbmc_umap <- gen_diffbin1_errors(highd_data = data_pbmc, nldr_data = umap_pbmc) |>
  dplyr::mutate(method = "UMAP_5_min_dist_0.8")

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_umap_5_min_dist_0.8.rds")

###########

## For tsne
tsne_pbmc <- read_rds("data/pbmc3k/pbmc_tsne_5.rds")
names(tsne_pbmc) <- c("emb1", "emb2", "ID")

error_pbmc_tsne <- gen_diffbin1_errors(highd_data = data_pbmc, nldr_data = tsne_pbmc) |>
  dplyr::mutate(method = "tsne_5")

write_rds(error_pbmc_tsne, "data/pbmc3k/error_pbmc_tsne_5.rds")

###########

## For tsne
tsne_pbmc <- read_rds("data/pbmc3k/pbmc_tsne_30.rds")
names(tsne_pbmc) <- c("emb1", "emb2", "ID")

error_pbmc_tsne <- gen_diffbin1_errors(highd_data = data_pbmc, nldr_data = tsne_pbmc) |>
  dplyr::mutate(method = "tsne_30")

write_rds(error_pbmc_tsne, "data/pbmc3k/error_pbmc_tsne_30.rds")

###########

## For phate
phate_pbmc <- read_rds("data/pbmc3k/pbmc_phate_5.rds")
names(phate_pbmc) <- c("emb1", "emb2", "ID")

error_pbmc_phate <- gen_diffbin1_errors(highd_data = data_pbmc, nldr_data = phate_pbmc) |>
  dplyr::mutate(method = "phate_5")

write_rds(error_pbmc_phate, "data/pbmc3k/error_pbmc_phate_5.rds")

###########

## For trimap
trimap_pbmc <- read_rds("data/pbmc3k/pbmc_trimap_12_4_3.rds")
names(trimap_pbmc) <- c("emb1", "emb2", "ID")

error_pbmc_trimap <- gen_diffbin1_errors(highd_data = data_pbmc, nldr_data = trimap_pbmc) |>
  dplyr::mutate(method = "trimap_12_4_3")

write_rds(error_pbmc_trimap, "data/pbmc3k/error_pbmc_trimap_12_4_3.rds")

###########

## For pacmap
pacmap_pbmc <- read_rds("data/pbmc3k/pbmc_pacmap_30_random_0.9_5.rds")
names(pacmap_pbmc) <- c("emb1", "emb2", "ID")

error_pbmc_pacamp <- gen_diffbin1_errors(highd_data = data_pbmc, nldr_data = pacmap_pbmc) |>
  dplyr::mutate(method = "pacmap")

write_rds(error_pbmc_pacamp, "data/pbmc3k/error_pbmc_pacmap_30_random_0.9_5.rds")

###########

## For umap
umap_pbmc <- read_rds("data/pbmc3k/pbmc_umap_9_min_dist_0.5.rds")
names(umap_pbmc) <- c("emb1", "emb2", "ID")

error_pbmc_umap <- gen_diffbin1_errors(highd_data = data_pbmc, nldr_data = umap_pbmc) |>
  dplyr::mutate(method = "UMAP_9_min_dist_0.5")

write_rds(error_pbmc_umap, "data/pbmc3k/error_pbmc_umap_9_min_dist_0.5.rds")
