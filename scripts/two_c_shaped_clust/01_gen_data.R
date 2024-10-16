library(tibble)
library(dplyr)
library(langevitour)
library(readr)

set.seed(20240110)

# Function to gen a curvilinear cluster in 4D space with an offset
gen_curv1_3d <- function(n) {
  # Bias sampling of 'a' but keep it within the same range (-0.5 to 0)
  a <- 3 * pi * (stats::rbeta(n = n, shape1 = 1, shape2 = 10) * -0.9)
  x2 <- sin(a)
  x1 <- 2.0 * stats::runif(n = n)
  x3 <- sign(a) * (cos(a) - 1)

  df <- tibble(
    x1 = x1,
    x2 = x2,
    x3 = x3
  )

  return(df)
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
}

# Simulate some s_curve_noise

sample_size <- 1000
curve1 <- gen_curv1_3d(n = sample_size)
langevitour(curve1)

curve2 <- gen_curv2_3d(n = sample_size)
langevitour(curve2)


# Apply an offset to one of the clusters to create a distance between them
offset <- c(2, 2, 2)  # Adjust these values to set the desired distance
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

write_rds(df, here::here("data/two_c_shaped_clust/two_c_shaped_clust_data.rds"))
