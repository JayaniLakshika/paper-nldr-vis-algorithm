library(tibble)
library(dplyr)
library(langevitour)
library(readr)

set.seed(20240110)

gen_curv1_3d <- function(n = 100) {

  x1 <- runif(n, 0, 2)
  x2 <- runif(n, 0, 3)
  x3 <- -(x1^3 + x2) + runif(n, 0, 0.5)

  df <- tibble(
    x1 = x1,
    x2 = x2,
    x3 = x3
  )

  df

}

gen_curv2_3d <- function(n) {
  a <- 3 * pi * stats::runif(n = n, min = -0.5, max = 0)
  x1 <- sin(a)
  x2 <- 2.0 * stats::runif(n = n)
  x3 <- sign(a) * (cos(a) - 1)

  df <- tibble(
    x1 = x1,
    x2 = x2,
    x3 = x3
  )

  df
}

# Simulate some s_curve_noise

sample_size <- 1000
curve2 <- gen_curv1_3d(n = sample_size)
langevitour(curve2)

curve1 <- gen_curv2_3d(n = sample_size)
langevitour(curve1)


# Apply an offset to one of the clusters to create a distance between them
offset <- c(1.5, 2.3, 1)  # Adjust these values to set the desired distance
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

write_rds(df, here::here("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_data.rds"))

