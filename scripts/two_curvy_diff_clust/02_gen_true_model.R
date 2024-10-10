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
true_model <- data.frame()

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
  true_model <- rbind(true_model, new_points)
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

names(true_model) <- c("x1", "x2", "x3", "ID")

# Visualize with langevitour
langevitour(true_model |> select(-ID),
            lineFrom = connections$from,
            lineTo = connections$to)


curvy2 <- true_model |>
  select(-ID) |>
  select(x2, x3, x1)

names(curvy2) <- paste0("x", 1:3)
# Apply an offset to one of the clusters to create a distance between them
offset <- c(2.5, 2.5, 2.5)  # Adjust these values to set the desired distance
curvy2 <- sweep(curvy2, 2, offset, "+") |>
  mutate(ID = NROW(true_model) + row_number())


connections_curvy2 <- connections
connections_curvy2$from <- connections$from + n_samples * num_y_points
connections_curvy2$to <- connections$to + n_samples * num_y_points


model_data <- bind_rows(true_model, curvy2)

model_data <- model_data |>
  dplyr::mutate(x4 = mean(runif(NROW(model_data), -0.05, 0.05)),
                x5 = mean(runif(NROW(model_data), -0.02, 0.02)),
                x6 = mean(runif(NROW(model_data), -0.1, 0.1)),
                x7 = mean(runif(NROW(model_data), -0.01, 0.01)))

connections_all <- bind_rows(connections, connections_curvy2)

# Visualize with langevitour
langevitour(model_data |> select(-ID),
            lineFrom = connections_all$from,
            lineTo = connections_all$to)

model_data <- model_data |>
  select(x1, x2, x3, x4, x5, x6, x7, ID)

write_rds(model_data, "data/two_curvy_diff_clust/two_curvy_diff_clust_true_model.rds")
write_rds(connections_all, "data/two_curvy_diff_clust/two_curvy_diff_clust_true_model_connections.rds")

