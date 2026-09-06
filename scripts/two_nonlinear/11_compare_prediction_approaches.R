### This script is to generate embeddings for new test data
library(dplyr)
library(tibble)
library(readr)
library(conflicted)

library(Rtsne)
library(umap) #predict for uwot not working
library(phateR)
library(reticulate)

set.seed(20240110)

conflicts_prefer(umap::umap)

use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

### Generate UMAP for new test data
d_2NC7_ts <- read_rds(here::here("data/two_nonlinear/test_two_non_linear_diff_shaped_close_clusters.rds"))

## tSNE (default)
perplexity <- 47

tSNE_fit <- d_2NC7_ts |>
  dplyr::select(where(is.numeric)) |>
  Rtsne::Rtsne(perplexity = perplexity,
               pca = FALSE)

## Loss function value
tail(tSNE_fit$itercosts, 1)

tSNE_data <- tSNE_fit$Y |>
  tibble::as_tibble(.name_repair = "unique")
names(tSNE_data) <- c("tSNE1", "tSNE2")

write_rds(tSNE_data, file = paste0("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_", perplexity, "_test.rds"))
