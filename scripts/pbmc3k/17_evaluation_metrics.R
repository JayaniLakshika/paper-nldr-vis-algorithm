library(tidyverse)
# Load the reticulate library
library(reticulate)
library(here)
library(FNN)
set.seed(20240110)
use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

# Run the Python script
py_run_file("~/Desktop/PhD Monash research files/Research papers/paper-nldr-vis-algorithm/scripts/evaluation.py")

# Create an instance of the Python class
my_instance <- py$MyClass()

## R_NX

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

highd_data <- read_rds("data/pbmc3k/pbmc_pca_50.rds")
names(highd_data) <- paste0("x", 1:50)

highd_data <- highd_data[, 1:9]

embeddings <- list(
  UMAP_5_min_dist_0.8     = read_rds(here("data/pbmc3k/pbmc_umap_5_min_dist_0.8.rds")),
  UMAP_30_min_dist_0.3     = read_rds(here("data/pbmc3k/pbmc_umap_30_min_dist_0.3.rds")),
  UMAP_5_min_dist_0.01     = read_rds(here("data/pbmc3k/pbmc_umap_5_min_dist_0.01.rds")),
  tsne_5     = read_rds(here("data/pbmc3k/pbmc_tsne_5.rds")),
  tsne_30     = read_rds(here("data/pbmc3k/pbmc_tsne_30.rds")),
  phate_5     = read_rds(here("data/pbmc3k/pbmc_phate_5.rds")),
  trimap_12_4_3    = read_rds(here("data/pbmc3k/pbmc_trimap_12_4_3.rds")),
  pacmap    = read_rds(here("data/pbmc3k/pbmc_pacmap_30_random_0.9_5.rds"))
)

evaluate_embedding <- function(name, lowd, highd) {
  # Global score from Python
  global_score <- my_instance$global_score(as.matrix(highd), as.matrix(lowd))

  # Random triplet accuracy (Python)
  triplet_acc <- py$evaluate_output(
    X = as.matrix(highd),
    X_new = as.matrix(lowd),
    y = NULL,
    name = name,
    baseline = FALSE,
    labelled = FALSE
  )$rte

  # R_NX AUC
  RNX_df <- R_NX(highd, lowd, max_k = 2620)
  RNX_AUC <- R_NX_AUC(RNX_df)

  # Shepard correlation (Spearman)
  shepard_df <- Shepard_diagram(highd, lowd)
  shepard_corr <- cor(shepard_df$dX, shepard_df$dY, method = "spearman")

  tibble(
    Method = name,
    Global_Score = global_score,
    Random_Triplet_Accuracy = triplet_acc,
    R_NX_AUC = RNX_AUC,
    Shepard_Correlation = shepard_corr
  )
}

results_all <- purrr::map2_dfr(names(embeddings), embeddings, ~evaluate_embedding(.x, .y, highd_data))
results_all <- results_all |>
  rename(c("GS" = "Global_Score",
           "RTA" = "Random_Triplet_Accuracy",
           "AUC(R_NX)" = "R_NX_AUC",
           "r" = "Shepard_Correlation"))

results_all <- results_all |>
  mutate(RGS = 1 - GS) #to reverse GS

results_all

write_rds(results_all, "data/pbmc3k/nldr_eval_metrics.rds")


RNX_curves <- purrr::imap_dfr(embeddings, function(lowd, name) {
  R_NX(as.matrix(highd_data), as.matrix(lowd), max_k = 2620) |> mutate(method = name)
}) |>
  mutate(method = factor(method, levels = c("UMAP_5_min_dist_0.8", "UMAP_30_min_dist_0.3", "UMAP_5_min_dist_0.01","tsne_5", "tsne_30", "phate_5", "trimap_12_4_3", "pacmap")))

write_rds(RNX_curves, "data/pbmc3k/RNX_curves.rds")

ggplot(RNX_curves, aes(x = K, y = R_NX, color = method)) +
  geom_line(linewidth = 1.2) +
  scale_x_log10() +
  scale_color_manual(
    values=c('#e41a1c','#ff7f00','#4daf4a',
             "#a65628",'#636363', '#984ea3')) +
  labs(y = expression(R[NX](k)), x = "Neighborhood size (k)",
       title = "Neighborhood Preservation (R_NX Curves)") +
  theme_minimal()

# Compute Shepard data for all embeddings
shepard_all <- purrr::map2_dfr(
  names(embeddings),
  embeddings,
  ~Shepard_diagram(as.matrix(highd_data), as.matrix(.y)) |>
    mutate(Method = .x)
)

write_rds(shepard_all, "data/pbmc3k/shepard_all.rds")


ggplot(shepard_all, aes(x = dX, y = dY)) +
  geom_point(alpha = 0.25, size = 0.6) +
  #geom_smooth(method = "loess", color = "red", se = FALSE) +
  facet_wrap(~ Method, scales = "free") +
  labs(
    x = "High-dimensional distances",
    y = "Low-dimensional distances",
    title = "Shepard Diagrams for Different NLDR Methods"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

library(GGally)

ggparcoord(
  data = results_all,
  columns = 2:5,           # numeric columns to plot
  groupColumn = 1,         # color by method
  scale = "uniminmax"      # normalize between 0 and 1 for fair comparison
) +
  scale_color_manual(
    values=c('#e41a1c','#ff7f00','#4daf4a',
             "#a65628",'#636363', '#984ea3')) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Parallel Coordinate Plot of NLDR Evaluation Metrics",
    y = "Scaled Metric Value",
    x = "Evaluation Metric"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

