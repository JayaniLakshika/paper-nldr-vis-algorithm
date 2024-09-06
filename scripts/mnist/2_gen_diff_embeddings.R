## Run only once
library(dplyr)
# library(snedata)
# library(ggflowchart)
# library(purrr) ## map function
# library(gridExtra) ## for grid.arrange
# library(rsample)
# # library(DT)
# library(ggbeeswarm)
library(ggplot2)
library(readr)
library(reticulate)

library(Rtsne)
library(umap)
library(phateR)
library(reticulate)
library(patchwork)

#library(grid)

set.seed(20240110)

use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

reticulate::source_python(here::here("scripts/function_scripts/Fit_PacMAP_code.py"))
reticulate::source_python(here::here("scripts/function_scripts/Fit_TriMAP_code.py"))

source(here::here("scripts/nldr_code.R"))

mnist_10_pcs_of_digit_1 <- read_rds("data/mnist/mnist_10_pcs_of_digit_1.rds")

# ## PaCMAP
#
# tem_dir <- tempdir()
#
# Fit_PacMAP_data(mnist_10_pcs_of_digit_1, tem_dir)
#
# path <- file.path(tem_dir, "mnist_10_pcs_of_digit_1_without_class.csv")
# path2 <- file.path(tem_dir, "dataset_3_PaCMAP_values.csv")
#
# Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.9, as.integer(2), path, path2)
#
# PacMAP_data <- read_csv(path2)
# PacMAP_data <- PacMAP_data |>
#   mutate(ID = 1:NROW(PacMAP_data))
#
# write_rds(PacMAP_data, file = "data/mnist/mnist_pacmap.rds")


##############
### tSNE
#tSNE_data_gau <- Fit_tSNE(mnist_10_pcs_of_digit_1, opt_perplexity = calculate_effective_perplexity(mnist_10_pcs_of_digit_1), with_seed = 20240110)
tSNE_data_gau <- Fit_tSNE(mnist_10_pcs_of_digit_1, opt_perplexity = 30, with_seed = 20240110)

tSNE_fit <- mnist_10_pcs_of_digit_1 |>
  Rtsne::Rtsne(perplexity = 30)

tSNE_data_gau <- tSNE_data_gau |>
  mutate(ID = 1:NROW(tSNE_data_gau))

plot_tSNE_2D(tSNE_data_gau)

#write_rds(tSNE_data_gau, file = paste0("data/mnist/mnist_tsne", calculate_effective_perplexity(mnist_10_pcs_of_digit_1),".rds"))
#write_rds(tSNE_data_gau, file = paste0("data/five_gau_clusters/tsne_data_five_gau_30.rds"))
write_rds(tSNE_data_gau, file = paste0("data/mnist/mnist_tsne30.rds"))


### UMAP

UMAP_fit <- umap(mnist_10_pcs_of_digit_1, n_neighbors = 15, n_components =  2)

UMAP_data_gau <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data_gau)[1:(ncol(UMAP_data_gau))] <- paste0(rep("UMAP",(ncol(UMAP_data_gau))), 1:(ncol(UMAP_data_gau)))

UMAP_data_gau <- UMAP_data_gau |>
  mutate(ID = 1:NROW(UMAP_data_gau))

write_rds(UMAP_data_gau, file = "data/mnist/mnist_umap.rds")

### TriMAP

# tem_dir <- tempdir()
#
# Fit_TriMAP_data(mnist_10_pcs_of_digit_1, tem_dir)
#
# path <- file.path(tem_dir, "mnist_10_pcs_of_digit_1_without_class.csv")
# path2 <- file.path(tem_dir, "dataset_3_TriMAP_values.csv")
#
# #Fit_TriMAP(as.integer(2), as.integer(5), as.integer(4), as.integer(3), path, path2)
# Fit_TriMAP(as.integer(2), as.integer(12), as.integer(4), as.integer(3), path, path2)
#
# TriMAP_data <- read_csv(path2)

trimap <- reticulate::import("trimap")

data_vector <- unlist(mnist_10_pcs_of_digit_1)
# Convert the vector into a matrix
data_matrix <- matrix(data_vector, ncol = NCOL(mnist_10_pcs_of_digit_1))

# Initialize PaCMAP instance
reducer <- trimap$TRIMAP(n_dims = as.integer(2),
                         n_inliers = as.integer(12),
                         n_outliers = as.integer(4),
                         n_random = as.integer(3))

# Perform dimensionality Reduction
TriMAP_data <- reducer$fit_transform(data_matrix) |>
  as_tibble()

names(TriMAP_data) <- c("TriMAP1", "TriMAP2")

TriMAP_data <- TriMAP_data |>
  mutate(ID = 1:NROW(TriMAP_data))
write_rds(TriMAP_data, file = "data/mnist/mnist_trimap.rds")


### PaCMAP

# tem_dir <- tempdir()
#
# Fit_PacMAP_data(mnist_10_pcs_of_digit_1, tem_dir)
#
# path <- file.path(tem_dir, "mnist_10_pcs_of_digit_1_without_class.csv")
# path2 <- file.path(tem_dir, "dataset_3_PaCMAP_values.csv")
#
# # Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.9, as.integer(2), path, path2)
# Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.5, as.integer(2), path, path2)
#
# PacMAP_data <- read_csv(path2)

pacmap <- reticulate::import("pacmap")

data_vector <- unlist(mnist_10_pcs_of_digit_1)
# Convert the vector into a matrix
data_matrix <- matrix(data_vector, ncol = NCOL(mnist_10_pcs_of_digit_1))

# Initialize PaCMAP instance
reducer <- pacmap$PaCMAP(n_components = as.integer(2), n_neighbors = as.integer(10), MN_ratio = 0.5, FP_ratio = as.integer(2))


# Perform dimensionality Reduction
PacMAP_data <- reducer$fit_transform(data_matrix, init = "random") |>
  as_tibble()

names(PacMAP_data) <- c("PaCMAP1", "PaCMAP2")

PacMAP_data <- PacMAP_data |>
  mutate(ID = 1:NROW(PacMAP_data))

write_rds(PacMAP_data, file = "data/mnist/mnist_pacmap.rds")

### Phate

PHATE_data <- Fit_PHATE(mnist_10_pcs_of_digit_1, knn = 5, with_seed = 20240110)
PHATE_data <- PHATE_data |>
  select(PHATE1, PHATE2)
PHATE_data <- PHATE_data |>
  mutate(ID = 1:NROW(PHATE_data))

write_rds(PHATE_data, file = "data/mnist/mnist_phate.rds")


