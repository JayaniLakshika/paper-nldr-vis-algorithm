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

gau_clust <- tibble::tibble(x1=rnorm(sample_size, mean = 0, sd = 0.2),
                            x2=rnorm(sample_size, mean = 0, sd = 0.2),
                            x3=rnorm(sample_size, mean = 1, sd = 0.2))

# Apply an offset to one of the clusters to create a distance between them
offset <- c(2, 2, 2)  # Adjust these values to set the desired distance
gau_clust <- sweep(gau_clust, 2, offset, "+")

data <- bind_rows(curve1,
                  gau_clust)

langevitour(data)

data$x4 <- runif(NROW(data), -0.05, 0.05)
data$x5 <- runif(NROW(data), -0.02, 0.02)
data$x6 <- runif(NROW(data), -0.1, 0.1)
data$x7 <- runif(NROW(data), -0.01, 0.01)

write_rds(data, here::here("data/one_curvy_one_gau_clust/one_curvy_one_gau_clust_data.rds"))


