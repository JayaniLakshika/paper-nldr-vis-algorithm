library(tibble)
library(dplyr)
library(readr)

set.seed(20260806)  # different seed from training (20240110)

# ── Same data-generating functions as training ──────────────────────────────

gen_curv1_3d <- function(n = 100) {
  x1 <- runif(n, 0, 2)
  x2 <- runif(n, 0, 3)
  x3 <- -(x1^3 + x2)
  x4 <- runif(n, 0, 1)
  tibble(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
}

gen_curv2_3d <- function(n) {
  a  <- 3 * pi * stats::runif(n = n, min = -0.5, max = 0)
  x1 <- sin(a)
  x2 <- 2.0 * stats::runif(n = n)
  x3 <- sign(a) * (cos(a) - 1)
  x4 <- cos(a)
  tibble(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
}

# ── Generate fresh test samples ──────────────────────────────────────────────

sample_size <- 250   # adjust if you want a different test set size

curve1 <- gen_curv2_3d(n = sample_size)
curve2 <- gen_curv1_3d(n = sample_size)

# Apply the same offset used in training
offset <- c(1.5, 3.3, 2, 1.5)
curve2 <- sweep(curve2, 2, offset, "+")

df_test <- bind_rows(curve1, curve2)

# Add the same noise dimensions as training
df_test$x5 <- runif(NROW(df_test), -0.02, 0.02)
df_test$x6 <- runif(NROW(df_test), -0.1,  0.1)
df_test$x7 <- runif(NROW(df_test), -0.01, 0.01)

# Add cluster labels (optional but useful for evaluation)
df_test$cluster <- rep(c("cluster1", "cluster2"), each = sample_size)

write_rds(
  df_test,
  here::here("data/two_nonlinear/test_two_non_linear_diff_shaped_close_clusters.rds")
)
