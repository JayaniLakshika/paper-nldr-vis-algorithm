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


projection <- cbind(
  c(0.5898,-0.4620,-0.1362,0.2047,-0.3680,0.1032,-0.4818),
  c(0.1756,0.0257,0.3683,0.2400,-0.5518,-0.4573,0.5115))
projected <- as.matrix(scaled_s_curve_noise) %*% projection

projected_df <- projected |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::rename(c("proj1" = "...1",
                  "proj2" = "...2")) |>
  dplyr::mutate(ID = dplyr::row_number())

limits <- 1
rng <- range(projected)
projected <- projected/max(abs(rng))
colnames(projected) <- c("P1", "P2")
projected <- data.frame(projected)
obs_labels <- as.character(1:nrow(training_data_scurve))

axis_scale <- limits/6
axis_pos <- -2/3 * limits

adj <- function(x) axis_pos + x * axis_scale
axes <- data.frame(x1 = adj(0),
                   y1 = adj(0),
                   x2 = adj(projection[, 1]),
                   y2 = adj(projection[, 2]))

axis_labels <- colnames(training_data_scurve |> select(-ID))
rownames(axes) <- axis_labels

theta <- seq(0, 2 * pi, length = 50)
circle <- data.frame(c1 = adj(cos(theta)), c2 = adj(sin(theta)))

projected_df |>
  ggplot(
    aes(
      x = proj1,
      y = proj2)) +
  geom_point(
    size = 0.8,
    alpha = 0.5,
    color = "#000000") +
  geom_segment(
    data=axes,
    aes(x=x1, y=y1, xend=x2, yend=y2),
    colour="grey70") +
  geom_text(
    data=axes,
    aes(x=x2, y=y2, label=rownames(axes)),
    colour="grey50",
    size = 5) +
  geom_path(
    data=circle,
    aes(x=c1, y=c2), colour="grey70")
