# Center the data by subtracting the mean of each column
center_data <- function(data) {
  apply(data, 2, function(col) col - mean(col))
}

# Function to scale data manually
scale_data_manual <- function(data) {
  # Step 1: Center the data (mean 0)
  data_centered <- center_data(data)

  # Step 2: Calculate the range of each dimension
  ranges <- apply(data_centered, 2, function(x) max(x) - min(x))

  # Step 3: Find the dimension with the largest range
  max_range <- max(ranges)
  scaled_data <- data_centered

  # Step 4: Scale each dimension
  for (i in seq_along(ranges)) {
    scaling_factor <- ranges[i] / max_range
    scaled_data[, i] <- data_centered[, i] / ranges[i] * scaling_factor
  }

  # Step 5: Rescale the dimension with the largest range to [-1, 1]
  largest_dim_index <- which.max(ranges)
  scaled_data[, largest_dim_index] <- scaled_data[, largest_dim_index] * max_range

  return(scaled_data)
}

training_data_scurve <- read_rds("data/s_curve/s_curve_training.rds")

# Apply the scaling function to the s_curve_noise data
scaled_s_curve_noise <- scale_data_manual(training_data_scurve |> select(-ID)) |>
  as_tibble()

langevitour(scaled_s_curve_noise, scale = 1)
