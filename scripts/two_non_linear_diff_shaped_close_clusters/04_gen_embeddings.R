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

#reticulate::source_python(paste0(here::here("R/function_scripts/Fit_PacMAP_code.py")))
#reticulate::source_python(paste0(here::here("R/function_scripts/Fit_TriMAP_code.py")))

#source(here::here("R/nldr_code.R"), local = TRUE)

data <- read_rds(here::here("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_data.rds"))

## Scale the data

# data <- data |>
#   mutate(across(everything(), ~ (. - mean(.)) / sd(.)))

## tSNE (default)
perplexity <- 30

tSNE_fit <- data |>
  dplyr::select(where(is.numeric)) |>
  Rtsne::Rtsne(perplexity = perplexity,
               pca = FALSE)

tSNE_data <- tSNE_fit$Y |>
  tibble::as_tibble(.name_repair = "unique")
names(tSNE_data) <- c("tSNE1", "tSNE2")

write_rds(tSNE_data, file = paste0("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_", perplexity, ".rds"))

## tSNE (another choice)
perplexity <- 62

tSNE_fit <- data |>
  dplyr::select(where(is.numeric)) |>
  Rtsne::Rtsne(perplexity = perplexity,
               pca = FALSE)

tSNE_data <- tSNE_fit$Y |>
  tibble::as_tibble(.name_repair = "unique")
names(tSNE_data) <- c("tSNE1", "tSNE2")

write_rds(tSNE_data, file = paste0("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_", perplexity, ".rds"))


## UMAP

n_neighbors <- 15
min_dist <- 0.1

# UMAP_model <- umap(data,
#                   n_neighbors = n_neighbors,
#                   min_dist = min_dist,
#                   n_components =  2,
#                   init ="spca")
#
# UMAP_data <- UMAP_model |>
#   as_tibble()

# UMAP_data <- UMAP_fit$layout |>
#   tibble::as_tibble(.name_repair = "unique")

# Create a config list with the desired parameters
umap_config <- umap.defaults
umap_config$n_neighbors <- n_neighbors      # Set the number of neighbors
umap_config$n_components <- 2    # Set the number of output dimensions (typically 2 or 3)
umap_config$min_dist <- min_dist

UMAP_fit <- umap(data, config = umap_config)

UMAP_data <- UMAP_fit$layout |>
  as_tibble()

names(UMAP_data) <- c("UMAP1", "UMAP2")

## Run only once
write_rds(UMAP_data, file = paste0("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_umap_n-neigbors_", n_neighbors, "_min-dist_", min_dist, ".rds"))

## Predict for true model
true_model_data <- read_rds("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_true_model.rds")
true_model_data <- true_model_data |>
  select(-ID)

predict_UMAP_df <- predict(UMAP_fit, true_model_data) |>
  as_tibble()

names(predict_UMAP_df) <- c("UMAP1", "UMAP2")

## Run only once
write_rds(predict_UMAP_df, file = "data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_umap_predict_true.rds")

## PHATE
knn <- 5

PHATE_data <- phate(data, knn = knn)
PHATE_data <- as_tibble(PHATE_data$embedding)

names(PHATE_data) <- c("PHATE1", "PHATE2")

write_rds(PHATE_data, file = paste0("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_phate_knn_", knn, ".rds"))


## TriMAP

trimap <- reticulate::import("trimap")

data_vector <- unlist(data)
# Convert the vector into a matrix
data_matrix <- matrix(data_vector, ncol = NCOL(data))

n_inliers <- as.integer(12)
n_outliers <- as.integer(4)
n_random <- as.integer(3)

# Initialize PaCMAP instance
reducer <- trimap$TRIMAP(n_dims = as.integer(2),
                         n_inliers = n_inliers,
                         n_outliers = n_outliers,
                         n_random = n_random)

# Perform dimensionality Reduction
TriMAP_data <- reducer$fit_transform(data_matrix) |>
  as_tibble()

names(TriMAP_data) <- c("TriMAP1", "TriMAP2")

write_rds(TriMAP_data, file = paste0("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_trimap_n-inliers_", n_inliers, "_n-outliers_", n_outliers, "_n-random_", n_random, ".rds"))

## PaCMAP

pacmap <- reticulate::import("pacmap")

data_vector <- unlist(data)
# Convert the vector into a matrix
data_matrix <- matrix(data_vector, ncol = NCOL(data))

n_neighbors <- as.integer(10)
MN_ratio <- 0.5
FP_ratio <- as.integer(2)
init <- "random"

# Initialize PaCMAP instance
reducer <- pacmap$PaCMAP(n_components = as.integer(2),
                         n_neighbors = n_neighbors,
                         MN_ratio = MN_ratio,
                         FP_ratio = FP_ratio)


# Perform dimensionality Reduction
PacMAP_data <- reducer$fit_transform(data_matrix, init = init) |>
  as_tibble()

names(PacMAP_data) <- c("PaCMAP1", "PaCMAP2")

write_rds(PacMAP_data, file = paste0("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_pacmap_n-neighbors_", n_neighbors,"_init_", init, "_MN-ratio_", MN_ratio, "_FP-ratio_", FP_ratio, ".rds"))

