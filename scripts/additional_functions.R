# Add plot label

interior_annotation <- function(label, position = c(0.9, 0.9), cex = 1, col="grey70") {
  annotation_custom(grid::textGrob(label = label,
                                   x = unit(position[1], "npc"), y = unit(position[2], "npc"),
                                   gp = grid::gpar(cex = cex, col=col)))
}

# Scale high-d data

# Center the data by subtracting the mean of each column
center_data <- function(data) {
  apply(data, 2, function(col) col - mean(col))
}

# Function to scale data manually
scale_data_manual <- function(data, type_col) {
  # Step 1: Center the data (mean 0)
  data_centered <- center_data(data |> select(-all_of(type_col)))

  # Step 2: Calculate the standard deviation of each dimension
  sds <- apply(data_centered, 2, sd)

  # Step 3: Scale each dimension to have the range [0, 1]
  data_scaled <- apply(data_centered, 2, function(col) col / max(abs(col)))

  # Step 4: Scale dimensions according to their variation
  # The dimension with the highest standard deviation is scaled to [-1, 1]
  # Other dimensions are scaled to smaller ranges based on their standard deviations
  max_sd <- max(sds)

  # Normalize the standard deviations to get scaling factors
  scaling_factors <- sds / max_sd

  for (i in seq_along(scaling_factors)) {
    data_scaled[, i] <- data_scaled[, i] * scaling_factors[i]
  }

  # Combine the scaled data with the 'type' column and return as a tibble
  data_scaled <- as_tibble(data_scaled) %>%
    mutate(!!type_col := data[[type_col]])

  return(data_scaled)
}

# Plot MSE

plot_mse <- function(error_df) {

  ggplot(error_df,
         aes(x = a1,
             y = RMSE,
             colour = method)) +
    geom_point(size = 0.8) +
    geom_line(linewidth = 0.3) +
    ylab("RMSE") +
    xlab(expression(paste("binwidth (", a[1], ")"))) +
    theme_minimal() +
    theme(panel.border = element_rect(fill = 'transparent'),
          plot.title = element_text(size = 12, hjust = 0.5, vjust = -0.5),
          axis.ticks.x = element_line(),
          axis.ticks.y = element_line(),
          legend.position = "none",
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          plot.margin = margin(0, 0, 0, 0))

}

# creating Standardization function
standardize = function(x){
  z <- (x - mean(x)) / sd(x)
  return( z)
}

# Plot digits in MNIST

plot_digit_img <- function(digit_df, palette, title_text, ncol = 5) {

  # Add a label column only to the first instance
  digit_df <- digit_df |>
    mutate(label = ifelse(instance == unique(instance)[1], title_text, NA))

  ggplot(data = digit_df, aes(x, y, fill = value)) +
    geom_tile() +
    geom_text(data = digit_df |> filter(!is.na(label)),
              aes(label = label),
              x = min(digit_df$x), y = max(digit_df$y), inherit.aes = FALSE,
              hjust = -0.3, vjust = 1.3, size = 7, color = "grey70") +
    facet_wrap(~ instance, ncol = ncol) +
    coord_fixed() +
    scale_fill_continuous_sequential(palette = palette) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0))

}

# Solve quadratic function

quad <- function(a = 3, b = 2 * a2, c = -(a2^2 + a1^2))
{
  a <- as.complex(a)
  answer <- c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
              (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
  if(all(Im(answer) == 0)) answer <- Re(answer)
  if(answer[1] == answer[2]) return(answer[1])
  answer[answer>0] ## only positive
}

plot_proj <- function(proj_obj,
                      point_param = c(1.5, 0.5, "#000000"), # size, alpha, color
                      line_param = c(0.5, 0.5, "#000000"), #linewidth, alpha
                      plot_limits,
                      axis_text_size = 3,
                      is_category = FALSE) {

  projected_df <- proj_obj$projected_df
  model_df <- proj_obj$model_df
  axes <- proj_obj$axes
  circle <- proj_obj$circle

  if(is_category == FALSE) {

    initial_plot <- projected_df |>
      ggplot(
        aes(
          x = proj1,
          y = proj2)) +
      geom_point(
        size = as.numeric(point_param[1]),
        alpha = as.numeric(point_param[2]),
        color = point_param[3]) +
      geom_segment(
        data = model_df,
        aes(
          x = proj1_from,
          y = proj2_from,
          xend = proj1_to,
          yend = proj2_to),
        color = line_param[3],
        linewidth = as.numeric(line_param[1]),
        alpha = as.numeric(line_param[2]))

  } else {

    projected_df <- projected_df |>
      dplyr::mutate(cluster = proj_obj$cluster)

    initial_plot <- projected_df |>
      ggplot(
        aes(
          x = proj1,
          y = proj2,
          colour = cluster)) +
      geom_point(
        size = as.numeric(point_param[1]),
        alpha = as.numeric(point_param[2])) +
      geom_segment(
        data = model_df,
        aes(
          x = proj1_from,
          y = proj2_from,
          xend = proj1_to,
          yend = proj2_to),
        color = line_param[3],
        linewidth = as.numeric(line_param[1]),
        alpha = as.numeric(line_param[2]))

  }

  initial_plot <- initial_plot +
    geom_segment(
      data=axes,
      aes(x=x1, y=y1, xend=x2, yend=y2),
      colour="grey70") +
    geom_text(
      data=axes,
      aes(x=x2, y=y2),
      label=rownames(axes),
      colour="grey50",
      size = axis_text_size) +
    geom_path(
      data=circle,
      aes(x=c1, y=c2), colour="grey70") +
    xlim(plot_limits) +
    ylim(plot_limits)

  initial_plot

}

