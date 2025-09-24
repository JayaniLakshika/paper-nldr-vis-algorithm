library(tidyverse)
# Load the reticulate library
library(reticulate)
use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

# Run the Python script
py_run_file("~/Desktop/PhD Monash research files/Research papers/paper-nldr-vis-algorithm/scripts/evaluation.py")

data <- read_rds(here::here("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_data.rds"))
tSNE_data <- read_rds(here::here("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_47.rds"))
# Create an instance of the Python class
my_instance <- py$MyClass()

score_result <- my_instance$global_score(as.matrix(data), as.matrix(tSNE_data))
score_result
results <- py$evaluate_output(X = as.matrix(data), X_new = as.matrix(tSNE_data), y = sample(0:4, 2000, replace = TRUE), name = "my_run", baseline = FALSE, labelled = TRUE)
results

## R_NX
library(FNN)      # for get.knn
library(dplyr)    # for pipes
library(ggplot2)  # for plotting

# X = high-d data (n x p)
# Y = low-d embedding (n x q)
R_NX_curve <- function(X, Y, k_max = 30) {
  n <- nrow(X)
  ks <- 1:k_max
  R_vals <- numeric(length(ks))

  # Precompute NN in high-D and low-D
  nn_X <- get.knn(X, k = k_max)$nn.index
  nn_Y <- get.knn(Y, k = k_max)$nn.index

  for (k in ks) {
    # Neighborhoods up to size k
    NX <- nn_X[, 1:k, drop = FALSE]
    NY <- nn_Y[, 1:k, drop = FALSE]

    # Compute preservation per point
    scores <- sapply(1:n, function(i) {
      length(intersect(NX[i, ], NY[i, ])) / k
    })

    R_vals[k] <- mean(scores)
  }

  tibble(k = ks, R_NX = R_vals)
}

# Example
set.seed(1)

curve_df <- R_NX_curve(as.matrix(data), as.matrix(tSNE_data), k_max = 20)

ggplot(curve_df, aes(x = k, y = R_NX)) +
  geom_line(linewidth=1.2, color="steelblue") +
  labs(y = expression(R[NX](k)), x = "Neighborhood size k") +
  theme_minimal()

## Shepard Diagram

Shepard_diagram <- function(X, Y, sample_pairs = 5000) {
  n <- nrow(X)

  # Compute all pairwise distances
  dX <- as.vector(dist(X))
  dY <- as.vector(dist(Y))

  # To avoid plotting all ~n^2 points, subsample
  if (length(dX) > sample_pairs) {
    idx <- sample(seq_along(dX), sample_pairs)
    dX <- dX[idx]
    dY <- dY[idx]
  }

  tibble(dX = dX, dY = dY)
}

# Example
shepard_df <- Shepard_diagram(as.matrix(data), as.matrix(tSNE_data))

ggplot(shepard_df, aes(x = dX, y = dY)) +
  geom_point(alpha=0.3, size=1) +
  geom_smooth(method="loess", color="red") +
  labs(x = "High-D distances", y = "Low-D distances",
       title = "Shepard Diagram") +
  theme_minimal()

cor(shepard_df$dX, shepard_df$dY, method = "spearman")

