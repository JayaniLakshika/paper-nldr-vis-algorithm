## This script is to check PCA results for Five Gaussian clusters data
library(tidyverse)
library(quollr)
library(patchwork)
library(cardinalR)

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

pca_ref_calc <- calculate_pca(data |> dplyr::select(where(is.numeric)))
data_pca <- pca_ref_calc$pca_components |>
  dplyr::mutate(ID = row_number()) |>
  dplyr::mutate(cluster = data$cluster)

rotations_df <- pca_ref_calc$rotations

