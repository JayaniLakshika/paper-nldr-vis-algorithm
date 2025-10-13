library(tidyverse)
# Load the reticulate library
library(reticulate)
set.seed(20240110)
use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

# Run the Python script
py_run_file("~/Desktop/PhD Monash research files/Research papers/paper-nldr-vis-algorithm/scripts/evaluation.py")

data <- read_rds(here::here("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_data.rds"))
tsne_two_curvy1 <- read_rds(here::here("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_47.rds"))
# Create an instance of the Python class
my_instance <- py$MyClass()

score_result <- my_instance$global_score(as.matrix(data), as.matrix(tsne_two_curvy1))
score_result
results <- py$evaluate_output(X = as.matrix(data), X_new = as.matrix(tsne_two_curvy1), y = NULL, name = "my_run", baseline = FALSE, labelled = FALSE)
results

## R_NX
library(FNN)

R_NX <- function(highd, lowd, max_k = 50) {
  n <- nrow(highd)
  Ks <- 1:max_k
  R_values <- numeric(length(Ks))

  # Get high-D and low-D neighbor indices
  nn_high <- get.knn(highd, k = max_k)$nn.index
  nn_low  <- get.knn(lowd, k = max_k)$nn.index

  for (k in Ks) {
    overlaps <- sapply(1:n, function(i) {
      length(intersect(nn_high[i, 1:k], nn_low[i, 1:k]))
    })

    Q_NX <- mean(overlaps / k)
    R_values[k] <- ((n - 1) * Q_NX - k) / (n - 1 - k)
  }

  tibble::tibble(K = Ks, R_NX = R_values)
}

R_NX_AUC <- function(res) {
  Ks <- res$K
  R <- res$R_NX
  sum(R / Ks) / sum(1 / Ks)
}

# Example
curve_df <- R_NX(as.matrix(data), as.matrix(tsne_two_curvy1), max_k = 20)

AUC_RNX <- R_NX_AUC(curve_df)
AUC_RNX

ggplot(curve_df, aes(x = K, y = R_NX)) +
  geom_line(linewidth=1.2, color="steelblue") +
  labs(y = expression(R[NX](k)), x = "Neighborhood size k") +
  theme_minimal()

## Shepard Diagram
## There is a R package: flipDimensionReduction and function GoodnessOfFitPlot()

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
shepard_df <- Shepard_diagram(as.matrix(data), as.matrix(tsne_two_curvy1))

ggplot(shepard_df, aes(x = dX, y = dY)) +
  geom_point(alpha=0.3, size=1) +
  geom_smooth(method="loess", color="red") +
  labs(x = "High-D distances", y = "Low-D distances",
       title = "Shepard Diagram") +
  theme_minimal()

cor(shepard_df$dX, shepard_df$dY, method = "spearman")

tsne_two_curvy2 <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_62.rds")
umap_two_curvy <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_umap_n-neigbors_15_min-dist_0.1.rds")
phate_two_curvy <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_phate_knn_5.rds")
trimap_two_curvy <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_trimap_n-inliers_12_n-outliers_4_n-random_3.rds")
pacmap_two_curvy <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_pacmap_n-neighbors_10_init_random_MN-ratio_0.5_FP-ratio_2.rds")

