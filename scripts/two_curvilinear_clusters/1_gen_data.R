# Function to generate a 2D curvilinear cluster with noise in the other dimensions
generate_curvilinear_2d_with_noise <- function(n, noise_dim1, noise_dim2) {
  x <- runif(n, 0, 2)
  y <- -(x^3 + runif(n, 0, 3)) + runif(n, 0, 0.5)

  # Add small noise to other dimensions
  z_noise <- rnorm(n, mean = 0, sd = noise_dim1)
  w_noise <- rnorm(n, mean = 0, sd = noise_dim2)

  return(cbind(x, y, z_noise, w_noise))
}

# Number of points for each cluster
n_points <- 500

# Generate the first curvilinear cluster in dimensions 1 and 2 with small noise in dimensions 3 and 4
curv_1_2 <- generate_curvilinear_2d_with_noise(n_points, noise_dim1 = 0.01, noise_dim2 = 0.01) |>
  as_tibble()

names(curv_1_2) <- c("x1", "x2", "x3", "x4")

# Generate the second curvilinear cluster in dimensions 3 and 4 with small noise in dimensions 1 and 2
curv_3_4 <- generate_curvilinear_2d_with_noise(n_points, noise_dim1 = 0.05, noise_dim2 = 0.05) |>
  as_tibble() |>
  select(z_noise, w_noise, x, y)

names(curv_3_4) <- c("x1", "x2", "x3", "x4")

# Apply an offset to one of the clusters to create a distance between them
offset <- c(2, 2, 2, 2)  # Adjust these values to set the desired distance
curv_3_4 <- sweep(curv_3_4, 2, offset, "+")

# Combine the clusters to create the final dataset
curvilinear_4d <- rbind(curv_1_2, curv_3_4)

# View the first few rows
head(curvilinear_4d)

# Visualize using langevitour
langevitour(curvilinear_4d)
