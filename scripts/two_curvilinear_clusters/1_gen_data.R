library(dplyr)

set.seed(20240110)

# Function to generate a 2D curvilinear cluster with noise in the other dimensions
generate_curvilinear_2d_with_noise <- function(n) {
  x1 <- runif(n, 0, 2)
  x2 <- -(x1^3 + runif(n, 0, 3)) + runif(n, 0, 0.5)

  # Add small noise to other dimensions
  x3 <- sin(x1 * pi) + runif(n, -0.5, 0.5)
  x4 <- rnorm(n, mean = 0, sd = 0.01)

  data <- tibble::tibble(x1 = x1,
                         x2 = x2,
                         x3 = x3,
                         x4 = x4)

  data
}

# Number of points for each cluster
n_points <- 500

# Generate the first curvilinear cluster in dimensions 1 and 2 with small noise in dimensions 3 and 4
curv_1_2 <- generate_curvilinear_2d_with_noise(n_points)

# Generate the second curvilinear cluster in dimensions 3 and 4 with small noise in dimensions 1 and 2
curv_3_4 <- generate_curvilinear_2d_with_noise(n_points) |>
  select(x3, x4, x1, x2)

names(curv_3_4) <- c("x1", "x2", "x3", "x4")

# Apply an offset to one of the clusters to create a distance between them
offset <- c(2, 2, 2, 2)  # Adjust these values to set the desired distance
curv_3_4 <- sweep(curv_3_4, 2, offset, "+")

# Combine the clusters to create the final dataset
curvilinear_4d <- bind_rows(curv_1_2, curv_3_4)

# View the first few rows
head(curvilinear_4d)

# Visualize using langevitour
langevitour(curvilinear_4d)

write_rds(curvilinear_4d, here::here("data/two_nonlinear_clusters/two_nonlinear_clusters_data.rds"))
