## Generate S-curve data with equal distant points

library(langevitour)
library(rsample)
library(scales)

set.seed(20230531)

## Assign n_samples and num_y_points 3:1
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

# Jitter parameters
#jitter_strength <- 0.01  # Adjust the jitter strength as needed

# Create a data frame for storing points
points <- data.frame()

# Loop through each point on the S-curve and generate thickness with jitter
for (i in 1:n_samples) {
  # Current (x, z) values
  current_x <- x[i]
  current_z <- z[i]

  # Generate y values with a thickness band
  y_values <- seq(0, 2, length.out = num_y_points)

  # Generate jitter only for y dimension
  #y_jitter <- rnorm(num_y_points, mean = 0, sd = jitter_strength)

  # Apply jitter to y dimension
  #y_values_with_jitter <- y_values + y_jitter

  # Create points for the current (x, z) with different y values
  new_points <- data.frame(x = rep(current_x, num_y_points),
                           y = y_values,
                           z = rep(current_z, num_y_points))

  # Combine with existing points
  points <- rbind(points, new_points)
}

# Visualize with langevitour
langevitour(points)

# Rename columns and add noise to other dimensions
names(points) <- c("x1", "x2", "x3")

# Sample size (same as the number of rows in the points data frame)
sample_size <- nrow(points)

# Add noise to other dimensions
points$x4 <- runif(sample_size, -0.02, 0.02)
points$x5 <- runif(sample_size, -0.02, 0.02)
points$x6 <- runif(sample_size, -0.01, 0.01)
points$x7 <- runif(sample_size, -0.01, 0.01)

points <- points |> mutate(ID = row_number())

write_rds(points, "data/s_curve/scurve_500.rds")

# # Add the ID to the s_curve_noise
# s_curve_noise <- points |>
#   dplyr::mutate(ID = dplyr::row_number()) |>
#   tibble::as_tibble()
#
# # Split the s_curve_noise as training and test
# data_split <- initial_split(s_curve_noise)
# s_curve_noise_training <- training(data_split) |>
#   dplyr::arrange(ID) |>
#   tibble::as_tibble()
#
# s_curve_noise_test <- testing(data_split) |>
#   dplyr::arrange(ID) |>
#   tibble::as_tibble()

#Combine with the true model for visualization
df <- dplyr::bind_rows(true_model, points |> mutate(type = "data"))

# Visualize with langevitour
langevitour(df |> dplyr::select(-type),
            lineFrom = append(point_connect_df$from, point_connect_df_dir_2$from),
            lineTo = append(point_connect_df$to, point_connect_df_dir_2$to),
            group = df$type,
            pointSize = append(rep(2, NROW(true_model)), rep(1, NROW(points))),
            levelColors = c("#6a3d9a", "#969696"),
            lineColors = rep("#969696", length(append(point_connect_df$to, point_connect_df_dir_2$from))))

