### This script is to generate an example for the effect of density

library(tibble)
library(dplyr)
library(langevitour)
library(readr)

set.seed(20240110)

# Function to generate C-shaped structure with diff density
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

sample_size <- 1000
curve1 <- gen_curv1_3d(n = sample_size)
langevitour(curve1)

write_rds(curve1, here::here("data/one_c_shaped_dens_structure/one_c_shaped_dens_data.rds"))

# Function to generate C-shaped structure with uniform density

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

sample_size <- 1000
curve2 <- gen_curv2_3d(n = sample_size)
langevitour(curve2)

write_rds(curve2, here::here("data/one_c_shaped_dens_structure/one_c_shaped_uni_dens_data.rds"))
