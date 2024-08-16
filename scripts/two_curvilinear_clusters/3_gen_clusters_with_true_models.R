library(dplyr)
library(tidyr)

# Function to generate a curvilinear pattern with noise in 2D and two noisy dimensions
generate_curvilinear_2d_with_noise <- function(n, noise_dim1, noise_dim2) {
  x <- runif(n, 0, 2)
  y <- -(x^3 + runif(n, 0, 3)) + runif(n, 0, 0.5)

  # Add small noise to other dimensions
  z_noise <- rnorm(n, mean = 0, sd = noise_dim1)
  w_noise <- rnorm(n, mean = 0, sd = noise_dim2)

  return(cbind(x, y, z_noise, w_noise))
}

# Function to generate a curvilinear grid pattern in 2D
generate_curvilinear_grid_2d <- function(n_grid_x, n_grid_y) {
  x <- seq(0, 2, length.out = n_grid_x)
  y <- seq(-0.5, 3, length.out = n_grid_y)

  # Use expand.grid to create a grid of (x, y) pairs
  curvilinear_grid <- expand.grid(x = x, y = y)

  # Adjust y values according to a curvilinear pattern based on x
  curvilinear_grid$y <- -(curvilinear_grid$x^3 + curvilinear_grid$y)

  return(curvilinear_grid)
}

# Generate the curvilinear data with noise
n_points <- 500
curvilinear_data1 <- generate_curvilinear_2d_with_noise(n_points, noise_dim1 = 0.01, noise_dim2 = 0.01) %>%
  as_tibble()

names(curvilinear_data1) <- c("x", "y", "z_noise", "w_noise")

curvilinear_data2 <- generate_curvilinear_2d_with_noise(n_points, noise_dim1 = 0.01, noise_dim2 = 0.01) %>%
  as_tibble()

names(curvilinear_data2) <- c("x", "y", "z_noise", "w_noise")

curvilinear_data2 <- curvilinear_data2 |>
  select("z_noise", "w_noise", "x", "y")

names(curvilinear_data2) <- c("x", "y", "z_noise", "w_noise")

# Apply an offset to one of the clusters to create a distance between them
offset <- c(2, 2, 2, 2)  # Adjust these values to set the desired distance
curvilinear_data2 <- sweep(curvilinear_data2, 2, offset, "+")

curvilinear_data <- bind_rows(curvilinear_data1, curvilinear_data2)

curvilinear_data <- curvilinear_data |>
  mutate(type = "data")

# Generate the grid data
n_grid_x <- 10  # Number of grid points along x-axis
n_grid_y <- 10  # Number of grid points along y-axis

curvilinear_grid1 <- generate_curvilinear_grid_2d(n_grid_x, n_grid_y) %>%
  as_tibble()

# Add small noise to the grid data to match dimensions with curvilinear data
curvilinear_grid1 <- curvilinear_grid1 %>%
  arrange(x, y) %>%
  mutate(
    z_noise = rnorm(n_grid_x * n_grid_y, mean = 0, sd = 0.01),
    w_noise = rnorm(n_grid_x * n_grid_y, mean = 0, sd = 0.01)
  ) |>
  mutate(ID = row_number())

curvilinear_grid2 <- generate_curvilinear_grid_2d(n_grid_x, n_grid_y) %>%
  as_tibble()

# Add small noise to the grid data to match dimensions with curvilinear data
curvilinear_grid2 <- curvilinear_grid2 %>%
  arrange(x, y) %>%
  mutate(
    z_noise = rnorm(n_grid_x * n_grid_y, mean = 0, sd = 0.01),
    w_noise = rnorm(n_grid_x * n_grid_y, mean = 0, sd = 0.01)
  ) |>
  select(z_noise, w_noise, x, y)

names(curvilinear_grid2) <- c("x", "y", "z_noise", "w_noise")

# Apply an offset to one of the clusters to create a distance between them
offset <- c(2, 2, 2, 2)  # Adjust these values to set the desired distance
curvilinear_grid2 <- sweep(curvilinear_grid2, 2, offset, "+") |>
  mutate(ID = 100 + row_number())

curvilinear_grid <- bind_rows(curvilinear_grid1,
                              curvilinear_grid2)

curvilinear_grid <- curvilinear_grid |>
  mutate(type = "model")

# Combine the curvilinear data and grid data
combined_data <- bind_rows(curvilinear_grid |> select(-ID), curvilinear_data)

# Create a dataframe for "from" and "to" connections within rows and columns
# Function to generate connections for a grid
generate_grid_connections <- function(grid_data, n_grid_x, n_grid_y) {
  # Generate connections between adjacent points
  connections <- data.frame()

  # Loop through rows
  for (i in 1:(n_grid_x * n_grid_y)) {
    # Calculate row and column indices
    row_idx <- ceiling(i / n_grid_x)
    col_idx <- i %% n_grid_x
    if (col_idx == 0) col_idx <- n_grid_x

    # Get the current point ID
    current_id <- grid_data$ID[i]

    # Connect to the right
    if (col_idx < n_grid_x) {
      right_id <- grid_data$ID[i + 1]
      connections <- rbind(connections, data.frame(From = current_id, To = right_id))
    }

    # Connect to the next row
    if (row_idx < n_grid_y) {
      below_id <- grid_data$ID[i + n_grid_x]
      connections <- rbind(connections, data.frame(From = current_id, To = below_id))
    }
  }

  return(connections)
}

# Generate connections for each grid
connections1 <- generate_grid_connections(curvilinear_grid1, n_grid_x, n_grid_y)
connections2 <- generate_grid_connections(curvilinear_grid2, n_grid_x, n_grid_y)

# Combine grid connections
all_connections <- bind_rows(connections1, connections2)

langevitour(combined_data |> select(-type),
            group = combined_data$type,
            lineFrom = all_connections$From,
            lineTo = all_connections$To)
