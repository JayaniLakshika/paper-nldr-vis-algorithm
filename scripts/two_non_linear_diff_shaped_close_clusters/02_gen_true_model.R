library(langevitour)
library(rsample)
library(scales)
library(dplyr)
library(readr)
library(tidyr)

set.seed(20230531)

## C-shaped curve
# Define the number of grid points per dimension
n_a <- 20  # Number of grid points for 'a'
n_x2 <- 5  # Number of grid points for 'x2'

# Create structured sequences for a and x2
a_seq <- seq(3 * pi * -0.5, 3 * pi * 0, length.out = n_a)
x2_seq <- seq(0, 2.0, length.out = n_x2)

# Generate a grid using expand.grid
grid_data1 <- expand.grid(a = a_seq, x2 = x2_seq)

# Compute the corresponding non-linear transformations
grid_data1$x1 <- sin(grid_data1$a)
grid_data1$x3 <- sign(grid_data1$a) * (cos(grid_data1$a) - 1)
grid_data1$x4 <- cos(grid_data1$a)

# View a subset of the grid
head(grid_data1)
langevitour::langevitour(grid_data1 |> dplyr::select(-a))

### Connections
# Add row index as an ID
grid_data1 <- grid_data1 |>
  mutate(ID = row_number())

grid_data1 <- grid_data1 |> dplyr::select(-a) |>
  select(x1, x2, x3, x4)

# grid_data1 <- grid_data1 |>
#   mutate(across(everything(), ~ (. - mean(.)) / sd(.)))

# Create connections between neighboring grid points
edges1 <- data.frame(from = integer(), to = integer())

# Define grid dimensions
n_a <- length(a_seq)
n_x2 <- length(x2_seq)

# Create connections
for (i in 1:(n_a * n_x2)) {
  row <- (i - 1) %% n_a + 1  # Row index
  col <- (i - 1) %/% n_a + 1  # Column index

  # Right neighbor
  if (row < n_a) {
    edges1 <- rbind(edges1, data.frame(from = i, to = i + 1))
  }

  # Top neighbor
  if (col < n_x2) {
    edges1 <- rbind(edges1, data.frame(from = i, to = i + n_a))
  }
}

# View some edges1
head(edges1)
langevitour::langevitour(grid_data1,
                         lineFrom = edges1$from,
                         lineTo = edges1$to)

## Non-linear curve
# Define grid resolution (number of points per dimension)
n_x1 <- 10  # Number of grid points for x1
n_x2 <- 10  # Number of grid points for x2
n_x4 <- 5   # Number of grid points for x4

# Generate evenly spaced sequences for x1, x2, x4
x1_seq <- seq(0, 2, length.out = n_x1)
x2_seq <- seq(0, 3, length.out = n_x2)
x4_seq <- seq(-0.5, 0.5, length.out = n_x4)

# Create a grid of (x1, x2, x4)
grid_data2 <- expand.grid(x1 = x1_seq, x2 = x2_seq, x4 = x4_seq)

# Compute x3 based on x1 and x2 (without random noise)
grid_data2$x3 <- -(grid_data2$x1^3 + grid_data2$x2)

# Add a unique row index (ID) to each grid point
grid_data2 <- grid_data2 |>
  mutate(ID = row_number())

grid_data2 <- grid_data2 |> dplyr::select(x1, x2, x3, x4)

# grid_data2 <- grid_data2 |>
#   mutate(across(everything(), ~ (. - mean(.)) / sd(.)))

# View a sample of the grid
head(grid_data2)

langevitour::langevitour(grid_data2)

### Connections

# Define grid dimensions
n_x1 <- length(x1_seq)
n_x2 <- length(x2_seq)
n_x4 <- length(x4_seq)

# Initialize empty data frame for edges2
edges2 <- data.frame(from = integer(), to = integer())

# Iterate through each point in the grid
for (i in 1:(n_x1 * n_x2 * n_x4)) {
  row <- (i - 1) %% n_x1 + 1     # x1 index
  col <- ((i - 1) %/% n_x1) %% n_x2 + 1  # x2 index
  depth <- (i - 1) %/% (n_x1 * n_x2) + 1  # x4 index

  # Right neighbor (x1 direction)
  if (row < n_x1) {
    edges2 <- rbind(edges2, data.frame(from = i, to = i + 1))
  }

  # Top neighbor (x2 direction)
  if (col < n_x2) {
    edges2 <- rbind(edges2, data.frame(from = i, to = i + n_x1))
  }

  # Front neighbor (x4 direction)
  if (depth < n_x4) {
    edges2 <- rbind(edges2, data.frame(from = i, to = i + (n_x1 * n_x2)))
  }
}

# View some edges2
head(edges2)

langevitour::langevitour(grid_data2,
                         lineFrom = edges2$from,
                         lineTo = edges2$to)

edges2 <- edges2 |>
  mutate(from = 100 + edges2$from,
         to = 100 + edges2$to)

# Apply an offset to one of the clusters to create a distance between them
offset <- c(1.5, 3.3, 2, 1.5)  # Adjust these values to set the desired distance
grid_data2 <- sweep(grid_data2, 2, offset, "+")

model_data <- bind_rows(grid_data1 |> dplyr::select(x1, x2, x3, x4),
                        grid_data2 |> dplyr::select(x1, x2, x3, x4))

connections_all <- bind_rows(edges1, edges2)

stats <- read_rds(here::here("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_data_stats.rds"))

model_data <- model_data |>
  mutate(across(everything(), ~ (. - stats[[paste0(cur_column(), "_mean")]]) /
                  stats[[paste0(cur_column(), "_sd")]]))

write_rds(model_data, "data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_true_model.rds")
write_rds(connections_all, "data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_true_model_connections.rds")

data <- read_rds(here::here("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_data_without_std.rds"))

df <- bind_rows(model_data |> mutate(type = "model"),
                data |> mutate(type = "data"))

# df <- df |>
#   mutate(across(-type, ~ (. - mean(.)) / sd(.)))

# Visualize with langevitour
langevitour(df |> select(-type),
            lineFrom = connections_all$from,
            lineTo = connections_all$to,
            group = df$type,
            levelColors = c("#6a3d9a", "#33a02c"))
