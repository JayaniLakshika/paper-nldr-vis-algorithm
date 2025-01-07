library(dplyr)
library(tibble)
library(readr)
library(conflicted)

library(Rtsne)
library(umap) #predit for uwot not working
library(phateR)
library(reticulate)

set.seed(20240110)

conflicts_prefer(umap::umap)

use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

### For uniform densed C-shaped structure

data <- read_rds(here::here("data/one_c_shaped_dens_structure/one_c_shaped_uni_dens_data.rds"))

## tSNE
perplexity <- 52

tSNE_fit <- data |>
  dplyr::select(where(is.numeric)) |>
  Rtsne::Rtsne(perplexity = perplexity,
               pca = FALSE)

tSNE_data <- tSNE_fit$Y |>
  tibble::as_tibble(.name_repair = "unique")
names(tSNE_data) <- c("tSNE1", "tSNE2")

write_rds(tSNE_data, file = paste0("data/one_c_shaped_dens_structure/one_c_shaped_uni_dens_structure_tsne_perplexity_", perplexity, ".rds"))

### For diff dens C-shaped structure

data <- read_rds(here::here("data/one_c_shaped_dens_structure/one_c_shaped_dens_data.rds"))

## tSNE
perplexity <- 52

tSNE_fit <- data |>
  dplyr::select(where(is.numeric)) |>
  Rtsne::Rtsne(perplexity = perplexity,
               pca = FALSE)

tSNE_data <- tSNE_fit$Y |>
  tibble::as_tibble(.name_repair = "unique")
names(tSNE_data) <- c("tSNE1", "tSNE2")

write_rds(tSNE_data, file = paste0("data/one_c_shaped_dens_structure/one_c_shaped_dens_structure_tsne_perplexity_", perplexity, ".rds"))
