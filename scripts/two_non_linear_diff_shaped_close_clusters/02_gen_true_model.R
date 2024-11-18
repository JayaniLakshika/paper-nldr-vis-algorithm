library(langevitour)
library(rsample)
library(scales)
library(dplyr)
library(readr)

set.seed(20230531)

# Number of points along the S-curve
n_samples <- 25

# Generate uniform parameter t
t <- seq(-0.5, 0, length.out = n_samples)

# Parametric equations for the S-curve in 3D
x <- sin(3 * pi * t)
z <- sign(t) * (cos(3 * pi * t) - 1)

# Set band thickness for the S-curve in the third dimension (y)
band_thickness <- 2.0
num_y_points <- 10  # Number of points along the y direction for each (x, z) pair

# Create a data frame for storing points
true_model1 <- data.frame()

# Loop through each point on the S-curve and generate thickness without jitter
for (i in 1:n_samples) {
  # Current (x, z) values
  current_x <- x[i]
  current_z <- z[i]

  # Generate y values with a thickness band
  y_values <- seq(0, band_thickness, length.out = num_y_points)

  # Create points for the current (x, z) with different y values
  new_points <- data.frame(x = rep(current_x, num_y_points),
                           y = y_values,
                           z = rep(current_z, num_y_points),
                           ID = seq((i-1)*num_y_points + 1, i*num_y_points))

  # Combine with existing points
  true_model1 <- rbind(true_model1, new_points)
}

# Create the connections (from-to pairs)
connections <- data.frame(
  from = integer(0),
  to = integer(0)
)

# Connect points along the length of the S-curve
for (i in 1:(n_samples - 1)) {
  for (j in 1:num_y_points) {
    from_id <- (i - 1) * num_y_points + j
    to_id <- i * num_y_points + j
    connections <- rbind(connections, data.frame(from = from_id, to = to_id))
  }
}

# Connect points along the width (y-direction)
for (i in 1:n_samples) {
  for (j in 1:(num_y_points - 1)) {
    from_id <- (i - 1) * num_y_points + j
    to_id <- from_id + 1
    connections <- rbind(connections, data.frame(from = from_id, to = to_id))
  }
}

names(true_model1) <- c("x1", "x2", "x3", "ID")

# Visualize with langevitour
langevitour(true_model1 |> select(-ID),
            lineFrom = connections$from,
            lineTo = connections$to)

# Function to generate a curvilinear grid pattern in 2D
generate_curvilinear_grid_2d <- function(n_grid_x, n_grid_y) {
  x <- seq(0, 2, length.out = n_grid_x)
  y <- seq(-3, 0.5, length.out = n_grid_y)

  # Use expand.grid to create a grid of (x, y) pairs
  curvilinear_grid <- expand.grid(x1 = x, x2 = y)

  # Adjust y values according to a curvilinear pattern based on x
  curvilinear_grid$x3 <- -(curvilinear_grid$x1^3 + curvilinear_grid$x2)

  return(curvilinear_grid)
}

# Generate the grid data
n_grid_x <- 20  # Number of grid points along x-axis
n_grid_y <- 20  # Number of grid points along y-axis


curvilinear_grid1 <- generate_curvilinear_grid_2d(n_grid_x, n_grid_y) |>
  as_tibble()

# Apply an offset to one of the clusters to create a distance between them
offset <- c(1, 1.5, 1)  # Adjust these values to set the desired distance
curvilinear_grid1 <- sweep(curvilinear_grid1, 2, offset, "+")

# Add small noise to the grid data to match dimensions with curvilinear data
true_model2 <- curvilinear_grid1 |>
  mutate(ID = row_number())

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
      connections <- rbind(connections, data.frame(from = current_id, to = right_id))
    }

    # Connect to the next row
    if (row_idx < n_grid_y) {
      below_id <- grid_data$ID[i + n_grid_x]
      connections <- rbind(connections, data.frame(from = current_id, to = below_id))
    }
  }

  return(connections)
}

# Generate connections for each grid
connections1 <- generate_grid_connections(true_model2, n_grid_x, n_grid_y) |>
  as_tibble()

connections1 <- connections1 |>
  mutate(from = 250 + connections1$from,
         to = 250 + connections1$to)

model_data <- bind_rows(true_model1, true_model2)

model_data <- model_data |>
  dplyr::mutate(x4 = mean(runif(NROW(model_data), -0.05, 0.05)),
                x5 = mean(runif(NROW(model_data), -0.02, 0.02)),
                x6 = mean(runif(NROW(model_data), -0.1, 0.1)),
                x7 = mean(runif(NROW(model_data), -0.01, 0.01)))

connections_all <- bind_rows(connections, connections1)

model_data <- model_data |>
  select(x1, x2, x3, x4, x5, x6, x7, ID)

# Visualize with langevitour
langevitour(model_data |> select(-ID),
            lineFrom = connections_all$from,
            lineTo = connections_all$to)

write_rds(model_data, "data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_true_model.rds")
write_rds(connections_all, "data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_true_model_connections.rds")


data <- read_rds(here::here("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_data.rds"))

df <- bind_rows(model_data |> select(-ID) |> mutate(type = "model"),
                data |> mutate(type = "data"))

# Visualize with langevitour
langevitour(df |> select(-type),
            lineFrom = connections_all$from,
            lineTo = connections_all$to,
            group = df$type,
            levelColors = c("#6a3d9a", "#33a02c"))
