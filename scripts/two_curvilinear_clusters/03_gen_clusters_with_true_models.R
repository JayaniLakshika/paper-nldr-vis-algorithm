library(dplyr)
library(tidyr)

# Function to generate a 2D curvilinear cluster with noise in the other dimensions
generate_curvilinear_2d_with_noise <- function(n) {
  x1 <- runif(n, -1, 2)
  x2 <- runif(n, -0.5, 3)
  x3 <- -(x1^3 + runif(n, 0, 3)) + runif(n, 0, 0.5)

  # Add small noise to other dimensions
  #x3 <- sin(x1 * pi) + runif(n, -0.5, 0.5)
  #x4 <- cos(x1 * pi) + runif(n, -0.5, 0.5)

  data <- tibble::tibble(x1 = x1,
                         x2 = x2,
                         x3 = x3) #,x4 = x4

  data
}

# Function to generate a curvilinear grid pattern in 2D
generate_curvilinear_grid_2d <- function(n_grid_x, n_grid_y) {
  x1 <- seq(-1, 2, length.out = n_grid_x)
  x2 <- seq(-0.5, 3, length.out = n_grid_y)

  # Use expand.grid to create a grid of (x, y) pairs
  curvilinear_grid <- expand.grid(x1 = x1, x2 = x2)

  # Adjust y values according to a curvilinear pattern based on x
  curvilinear_grid$x3 <- -(curvilinear_grid$x1^3 + curvilinear_grid$x1)

  return(curvilinear_grid)
}

# Generate the curvilinear data with noise
n_points <- 750

# Generate the first curvilinear cluster in dimensions 1 and 2 with small noise in dimensions 3 and 4
curv_1_2 <- generate_curvilinear_2d_with_noise(n_points)

# Generate the second curvilinear cluster in dimensions 3 and 4 with small noise in dimensions 1 and 2
curv_3_4 <- generate_curvilinear_2d_with_noise(n_points) |>
  select(x3, x1, x2)

names(curv_3_4) <- c("x1", "x2", "x3")

# Apply an offset to one of the clusters to create a distance between them
offset <- c(4, 4, 4)  # Adjust these values to set the desired distance
curv_3_4 <- sweep(curv_3_4, 2, offset, "+")

# Combine the clusters to create the final dataset
curvilinear_4d <- bind_rows(curv_1_2, curv_3_4)

curvilinear_data <- curvilinear_4d |>
  mutate(type = "data")

# Generate the grid data
n_grid_x <- 20  # Number of grid points along x-axis
n_grid_y <- 20  # Number of grid points along y-axis

curvilinear_grid1 <- generate_curvilinear_grid_2d(n_grid_x, n_grid_y) |>
  as_tibble()

# Add small noise to the grid data to match dimensions with curvilinear data
curvilinear_grid1 <- curvilinear_grid1 |>
  # arrange(x1, x2) |>
  # mutate(
  #   z_noise = rnorm(n_grid_x * n_grid_y, mean = 0, sd = 0.01),
  #   w_noise = rnorm(n_grid_x * n_grid_y, mean = 0, sd = 0.01)
  # ) |>
  mutate(ID = row_number())

curvilinear_grid2 <- generate_curvilinear_grid_2d(n_grid_x, n_grid_y) |>
  as_tibble()

# Add small noise to the grid data to match dimensions with curvilinear data
curvilinear_grid2 <- curvilinear_grid2 |>
  # arrange(x, y) |>
  # mutate(
  #   z_noise = rnorm(n_grid_x * n_grid_y, mean = 0, sd = 0.01),
  #   w_noise = rnorm(n_grid_x * n_grid_y, mean = 0, sd = 0.01)
  # ) |>
  select(x3, x1, x2)

names(curvilinear_grid2) <- c("x1", "x2", "x3")

# Apply an offset to one of the clusters to create a distance between them
offset <- c(4, 4, 4)  # Adjust these values to set the desired distance
curvilinear_grid2 <- sweep(curvilinear_grid2, 2, offset, "+") |>
  mutate(ID = n_grid_x * n_grid_y + row_number())

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
