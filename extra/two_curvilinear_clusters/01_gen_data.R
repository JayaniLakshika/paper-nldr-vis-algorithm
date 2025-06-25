library(tibble)
library(dplyr)
library(langevitour)
library(readr)

set.seed(20240110)

curve_3d <- function(n_samples = 100) {
  a <- 3 * pi * stats::runif(n = n_samples, min = -0.5, max = 0)
  x1 <- sin(a)
  x2 <- 2.0 * stats::runif(n = n_samples)
  x3 <- sign(a) * (cos(a) - 1)

  df <- tibble(
    x1 = x1,
    x2 = x2,
    x3 = x3
  )
}

# Simulate some s_curve_noise

sample_size <- 1000
curve1 <- curve_3d(n_samples = sample_size)
langevitour(curve1)

curve2 <- curve1 |>
  select(x2, x3, x1)

names(curve2) <- paste0("x", 1:3)

# Apply an offset to one of the clusters to create a distance between them
offset <- c(2.5, 2.5, 2.5)  # Adjust these values to set the desired distance
curve2 <- sweep(curve2, 2, offset, "+")

df <- bind_rows(
  curve1,
  curve2
)

df$x4 <- runif(NROW(df), -0.05, 0.05)
df$x5 <- runif(NROW(df), -0.02, 0.02)
df$x6 <- runif(NROW(df), -0.1, 0.1)
df$x7 <- runif(NROW(df), -0.01, 0.01)

langevitour(df)

write_rds(df, here::here("data/two_curvy_clust/two_curvy_clust_data.rds"))

