library(dplyr)
library(tidyr)

# Function to generate a curvilinear grid pattern in 2D
generate_curvilinear_grid_2d <- function(n_grid_x, n_grid_y) {
  x <- seq(0, 2, length.out = n_grid_x)
  y <- seq(-3, 0.5, length.out = n_grid_y)

  # Use expand.grid to create a grid of (x, y) pairs
  curvilinear_grid <- expand.grid(x = x, y = y)

  # Adjust y values according to a curvilinear pattern based on x
  curvilinear_grid$y <- -(curvilinear_grid$x^3 + curvilinear_grid$y)

  return(curvilinear_grid)
}

# Number of grid points in each dimension
n_grid_x <- 10  # Adjust this number for a finer or coarser grid along x
n_grid_y <- 10  # Adjust this number for a finer or coarser grid along y

# Generate the curvilinear grid of points in 2D
curvilinear_grid_2d <- generate_curvilinear_grid_2d(n_grid_x, n_grid_y)

# Add an ID column to uniquely identify each point
curvilinear_grid_2d <- curvilinear_grid_2d %>%
  arrange(x, y) %>%
  mutate(ID = row_number()) %>%
  mutate(
    z_noise = rnorm(n_grid_x * n_grid_y, mean = 0, sd = 0.01),
    w_noise = rnorm(n_grid_x * n_grid_y, mean = 0, sd = 0.02)
  )

# Create a dataframe for "from" and "to" connections within rows and columns
connections <- curvilinear_grid_2d %>%
  mutate(Next_ID_x = lead(ID, n = 1),
         Next_ID_y = lead(ID, n = n_grid_x)) %>%
  filter(!is.na(Next_ID_x) | !is.na(Next_ID_y)) %>%

  # Remove connections that wrap around to the next row or column
  mutate(Last_In_Row = (ID %% n_grid_x == 0),
         Last_In_Column = (ID > (n_grid_x * (n_grid_y - 1)))) %>%
  filter(!(Last_In_Row & !is.na(Next_ID_x)) &
           !(Last_In_Column & !is.na(Next_ID_y))) %>%

  select(From = ID, To_x = Next_ID_x, To_y = Next_ID_y) %>%
  pivot_longer(cols = c(To_x, To_y), names_to = "Direction", values_to = "To") %>%
  filter(!is.na(To)) %>%
  select(From, To)

# Identify the last grid points in each row and connect them sequentially
last_points <- curvilinear_grid_2d %>%
  filter(ID %% n_grid_x == 0) %>%
  arrange(ID)

# Create sequential connections between the last points in each row
last_point_connections <- tibble(
  From = last_points$ID[-nrow(last_points)],
  To = last_points$ID[-1]
)

# Combine all connections
all_connections <- bind_rows(connections, last_point_connections)

# View the combined connections dataframe
print(all_connections)

# Visualize the curvilinear grid and connections using langevitour
langevitour(
  curvilinear_grid_2d %>% select(-ID),
  lineFrom = all_connections$From,
  lineTo = all_connections$To
)

data <- bind_rows(curvilinear_grid_2d %>% select(-ID) |> mutate(type = "model"),
                  curv_1_2 |> mutate(type = "data"))
langevitour(
  data |> select(-type),
  group = data$type,
  lineFrom = all_connections$From,
  lineTo = all_connections$To
)
