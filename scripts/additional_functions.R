# Add plot label

interior_annotation <- function(label, position = c(0.92, 0.92), cex = 1, col="grey70") {
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
             y = sqrt(MSE),
             colour = method)) +
    geom_point(size = 0.8) +
    geom_line(linewidth = 0.3) +
    #scale_y_log10() +
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

plot_digit_img <- function(digit_df, palette, title_text, ncol = 3) {

  # Add a label column only to the first instance
  digit_df <- digit_df |>
    mutate(label = ifelse(instance == unique(instance)[1], title_text, NA))

  ggplot(data = digit_df, aes(x, y, fill = value)) +
    geom_tile() +
    geom_text(data = digit_df |> filter(!is.na(label)),
              aes(label = label),
              x = min(digit_df$x), y = max(digit_df$y), inherit.aes = FALSE,
              hjust = 0, vjust = 1, size = 10, color = "grey70") +
    facet_wrap(~ instance, ncol = ncol) +
    #ggtitle(title_text) +
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


