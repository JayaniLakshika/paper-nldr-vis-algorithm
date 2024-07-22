library(langevitour)
library(rsample)
library(scales)
library(dplyr)

set.seed(20230531)

# Number of points along the S-curve
n_samples <- 50

# Generate uniform parameter t
t <- seq(-0.5, 0.5, length.out = n_samples)

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

# Visualize with langevitour
langevitour(true_model)

# Rename columns and add noise to other dimensions
names(true_model) <- c("x1", "x2", "x3", "ID")

true_model <- true_model |>
  dplyr::mutate(x4 = mean(runif(NROW(true_model), -0.02, 0.02)),
                x5 = mean(runif(NROW(true_model), -0.02, 0.02)),
                x6 = mean(runif(NROW(true_model), -0.1, 0.1)),
                x7 = mean(runif(NROW(true_model), -0.01, 0.01)))

training_data <- read_rds("data/s_curve/s_curve_training.rds") |>
  dplyr::select(-ID) |>
  dplyr::mutate(type = "data")

# Combine with the true model for visualization
df <- dplyr::bind_rows(true_model |> select(-ID) |> mutate(type = "true model"), training_data)

# Visualize with langevitour
langevitour(df |> dplyr::select(-type),
            lineFrom = connections$from,
            lineTo = connections$to,
            group = df$type,
            pointSize = append(rep(2, NROW(true_model)), rep(1, NROW(training_data))),
            levelColors = c("#6a3d9a", "#969696"),
            lineColors = rep("#969696", nrow(connections)))
