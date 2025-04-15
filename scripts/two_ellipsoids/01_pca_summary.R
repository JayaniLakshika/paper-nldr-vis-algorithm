## This script is to check PCA results for Five Gaussian clusters data
library(tidyverse)
library(quollr)
library(patchwork)
library(cardinalR)

library(dplyr)
library(purrr) ## map function
library(rsample)
library(ggplot2)
library(readr)
library(geozoo)
library(mvtnorm)

library(Rtsne)
library(umap)
library(phateR)
library(reticulate)
library(patchwork)

library(grid)


set.seed(20240110)

use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

reticulate::source_python(here::here("scripts/function_scripts/Fit_PacMAP_code.py"))
reticulate::source_python(here::here("scripts/function_scripts/Fit_TriMAP_code.py"))

source(here::here("scripts/nldr_code.R"))

source("scripts/additional_functions.R")
set.seed(20240110)

clr_choice <- "#0077A3"

calculate_pca <- function(feature_dataset){
  pcaY_cal <- prcomp(feature_dataset, center = TRUE, scale = FALSE)
  PCAresults <- data.frame(pcaY_cal$x[, 1:4])
  pcaRotations <- data.frame(pcaY_cal$rotation[, 1:4])
  summary_pca <- summary(pcaY_cal)
  var_explained_df <- data.frame(PC= paste0("PC",1:4),
                                 var_explained=(pcaY_cal$sdev[1:2])^2/sum((pcaY_cal$sdev[1:2])^2))
  return(list(prcomp_out = pcaY_cal,pca_components = PCAresults, summary = summary_pca, var_explained_pca  = var_explained_df, rotations = pcaRotations))
}

# Stretch factors along each dimension (larger = more stretch)
stretch_factors <- c(2, 0.5, 2, 0.5)

# Create diagonal covariance matrix
Sigma <- diag(stretch_factors^2)

df1 <- gen_gaussian(n = 500, p = 4, m = c(0, 0, 0, 0), s = Sigma)
df2 <- gen_gaussian(n = 500, p = 4, m = c(0.5 * 6 ^2, 0, 0, 0), s = Sigma)

data <- dplyr::bind_rows(df1, df2)
#langevitour::langevitour(data)

### Generate embeddings (Run only once)

### tSNE
tSNE_data_gau <- Fit_tSNE(data, opt_perplexity = calculate_effective_perplexity(data), with_seed = 20240110)

tSNE_data_gau <- tSNE_data_gau |>
  select(-ID)

plot_tSNE_2D(tSNE_data_gau)

write_rds(tSNE_data_gau, file = paste0("data/two_ellipsoids/tsne_data_two_ellipsoids_", calculate_effective_perplexity(data),".rds"))
#write_rds(tSNE_data_gau, file = paste0("data/two_ellipsoids/tsne_data_two_ellipsoids_30.rds"))


### UMAP

UMAP_fit <- umap(data, n_neighbors = 15, n_components =  2)

UMAP_data_gau <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data_gau)[1:(ncol(UMAP_data_gau))] <- paste0(rep("UMAP",(ncol(UMAP_data_gau))), 1:(ncol(UMAP_data_gau)))

write_rds(UMAP_data_gau, file = "data/two_ellipsoids/umap_data_two_ellipsoids.rds")

### TriMAP

tem_dir <- tempdir()

Fit_TriMAP_data(data, tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_TriMAP_values.csv")

#Fit_TriMAP(as.integer(2), as.integer(5), as.integer(4), as.integer(3), path, path2)
Fit_TriMAP(as.integer(2), as.integer(12), as.integer(4), as.integer(3), path, path2)

TriMAP_data <- read_csv(path2)

write_rds(TriMAP_data, file = "data/two_ellipsoids/trimap_data_two_ellipsoids.rds")


### PaCMAP

tem_dir <- tempdir()

Fit_PacMAP_data(data, tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_PaCMAP_values.csv")

# Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.9, as.integer(2), path, path2)
Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.5, as.integer(2), path, path2)

PacMAP_data <- read_csv(path2)
write_rds(PacMAP_data, file = "data/two_ellipsoids/pacmap_data_two_ellipsoids.rds")



### Phate

PHATE_data <- Fit_PHATE(data, knn = 5, with_seed = 20240110)
PHATE_data <- PHATE_data |>
  select(PHATE1, PHATE2)

write_rds(PHATE_data, file = "data/two_ellipsoids/phate_data_two_ellipsoids.rds")



### Fit PCA
pca_ref_calc <- calculate_pca(data |> dplyr::select(where(is.numeric)))
data_pca <- pca_ref_calc$pca_components |>
  dplyr::mutate(ID = row_number()) |>
  dplyr::mutate(cluster = data$cluster)

rotations_df <- pca_ref_calc$rotations

