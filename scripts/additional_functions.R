# Add plot label

interior_annotation <- function(label, position = c(0.92, 0.92), cex = 1, col="grey70") {
  annotation_custom(grid::textGrob(label = label,
                                   x = unit(position[1], "npc"), y = unit(position[2], "npc"),
                                   gp = grid::gpar(cex = cex, col=col)))
}

# Generate the high-d-vis obj

gen_nldr_vis_algo_obj <- function(high_d_data, nldr_data, num_x_bins) {

  nldr_scaled_obj <- gen_scaled_data(
    data = nldr_data)

  nldr_scaled <- nldr_scaled_obj$scaled_nldr |>
    mutate(ID = row_number())

  lim1 <- nldr_scaled_obj$lim1
  lim2 <- nldr_scaled_obj$lim2
  r2 <- diff(lim2)/diff(lim1)

  ## hexagon binning to have regular hexagons
  hb_obj <- hex_binning(
    data = nldr_scaled,
    bin1 = num_x_bins,
    r2 = r2,
    q = 0.1)

  ## Data set with all centroids
  all_centroids_df <- hb_obj$centroids

  a1 <- hb_obj$a1

  ## Generate all coordinates of hexagons
  hex_grid <- hb_obj$hex_poly

  ## To obtain the standardise counts within hexbins
  counts_df <- hb_obj$std_cts
  df_bin_centroids <- extract_hexbin_centroids(
    centroids_df = all_centroids_df,
    counts_df = counts_df) |>
    filter(drop_empty == FALSE)

  nldr_with_hb_id <- hb_obj$data_hb_id
  df_all <- dplyr::bind_cols(high_d_data,
                             nldr_with_hb_id)
  df_bin <- avg_highd_data(data = df_all)

  hex_grid_with_counts <-
    left_join(hex_grid,
              counts_df,
              by = c("hex_poly_id" = "hb_id"))

  hex_grid_nonempty <- hex_grid |>
    filter(hex_poly_id %in% df_bin_centroids$hexID)

  tr1_object <- tri_bin_centroids(
    df_bin_centroids, x = "c_x", y = "c_y")
  tr_from_to_df <- gen_edges(
    tri_object = tr1_object)

  ## Compute 2D distances
  distance_df <- cal_2d_dist(
    tr_coord_df = tr_from_to_df,
    start_x = "x_from",
    start_y = "y_from",
    end_x = "x_to",
    end_y = "y_to",
    select_vars = c("from", "to", "distance"))

  ## To find the benchmark value
  benchmark <- find_lg_benchmark(
    distance_edges = distance_df,
    distance_col = "distance")

  return(list(nldr_scaled = nldr_scaled,
              a1 = a1,
              hex_grid_with_counts = hex_grid_with_counts,
              hex_grid_nonempty = hex_grid_nonempty,
              df_bin_centroids = df_bin_centroids,
              tr_from_to_df = tr_from_to_df,
              df_bin = df_bin,
              distance_df = distance_df,
              benchmark = benchmark))

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

# Get projection

get_projection <- function(projection, proj_scale, scaled_data,
                           scaled_data_model, distance_df_small_edges,
                           axis_param) {

  projection_scaled <- projection * proj_scale
  projected <- as.matrix(scaled_data) %*% projection_scaled

  projected_df <- projected |>
    tibble::as_tibble(.name_repair = "unique") |>
    dplyr::rename(c("proj1" = "...1",
                    "proj2" = "...2")) |>
    dplyr::mutate(ID = dplyr::row_number())

  projected_model <- as.matrix(scaled_data_model) %*% projection_scaled

  projected_model_df <- projected_model |>
    tibble::as_tibble(.name_repair = "unique") |>
    dplyr::rename(c("proj1" = "...1",
                    "proj2" = "...2")) |>
    dplyr::mutate(ID = dplyr::row_number())

  model_df <- dplyr::left_join(
    distance_df_small_edges |> select(-distance),
    projected_model_df,
    by = c("from" = "ID"))

  names(model_df)[3:NCOL(model_df)] <- paste0(names(projected_model_df)[-NCOL(projected_model_df)], "_from")

  model_df <- dplyr::left_join(model_df, projected_model_df, by = c("to" = "ID"))
  names(model_df)[(2 + NCOL(projected_model_df)):NCOL(model_df)] <- paste0(names(projected_model_df)[-NCOL(projected_model_df)], "_to")

  limits <- axis_param$limits
  axis_scaled <- axis_param$axis_scaled
  axis_pos_x <- axis_param$axis_pos_x
  axis_pos_y <- axis_param$axis_pos_y
  threshold <- axis_param$threshold

  axes_obj <- gen_axes(
    proj = projection * axis_scaled,
    limits = limits,
    axis_pos_x = axis_pos_x,
    axis_pos_y = axis_pos_y,
    axis_labels = names(scaled_data),
    threshold = threshold)

  axes <- axes_obj$axes
  circle <- axes_obj$circle

  return(list(projected_df = projected_df,
              model_df = model_df,
              axes = axes,
              circle = circle))

}

# Plot projection

plot_proj <- function(proj_obj,
                      point_param = c(1.5, 0.5, "#000000"), # size, alpha, color
                      line_param = c(0.5, 0.5, "#000000"), #linewidth, alpha
                      plot_limits, title, cex = 2,
                      position = c(0.92, 0.92),
                      axis_text_size = 3,
                      is_color = FALSE) {

  projected_df <- proj_obj$projected_df
  model_df <- proj_obj$model_df
  axes <- proj_obj$axes
  circle <- proj_obj$circle

  if(is_color == FALSE) {

    initial_plot <- ggplot() +
      geom_point(
        data = projected_df,
        aes(
          x = proj1,
          y = proj2),
        size = as.numeric(point_param[1]),
        alpha = as.numeric(point_param[2]),
        color = point_param[3])

  } else {

    projected_df <- projected_df |>
      dplyr::mutate(cluster = proj_obj$cluster)

    initial_plot <- ggplot() +
      geom_point(
        data = projected_df,
        aes(
          x = proj1,
          y = proj2,
          color = cluster),
        size = as.numeric(point_param[1]),
        alpha = as.numeric(point_param[2]))

  }

  initial_plot <- initial_plot +
    geom_segment(
      data = model_df,
      aes(
        x = proj1_from,
        y = proj2_from,
        xend = proj1_to,
        yend = proj2_to),
      color = line_param[3],
      linewidth = as.numeric(line_param[1]),
      alpha = as.numeric(line_param[2])) +
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
    ylim(plot_limits) +
    interior_annotation(title, position, cex = cex)

  initial_plot

}

# Generate axes

gen_axes <- function(proj, limits = 1, axis_pos_x = NULL, axis_pos_y = NULL, axis_labels, threshold) {

  axis_scale <- limits/6

  if (is.null(axis_pos_x)) {

    axis_pos_x <- -2/3 * limits

  }

  if (is.null(axis_pos_y)) {

    axis_pos_y <- -2/3 * limits

  }

  adj <- function(x, axis_pos) axis_pos + x * axis_scale
  axes <- data.frame(x1 = adj(0, axis_pos_x),
                     y1 = adj(0, axis_pos_y),
                     x2 = adj(proj[, 1], axis_pos_x),
                     y2 = adj(proj[, 2], axis_pos_y))

  rownames(axes) <- axis_labels

  ## To remove axes
  axes <- axes |>
    mutate(distance = sqrt((x2 - x1)^2 + (y2 - y1)^2)) |>
    filter(distance >= threshold)

  theta <- seq(0, 2 * pi, length = 50)
  circle <- data.frame(c1 = adj(cos(theta), axis_pos_x),
                       c2 = adj(sin(theta), axis_pos_y))

  return(list(axes = axes, circle = circle))

}

# Plot MSE

plot_mse <- function(error_df) {

  ggplot(error_df,
         aes(x = a1,
             y = MSE,
             colour = method)) +
    geom_point(size = 0.8) +
    geom_line(linewidth = 0.3) +
    scale_y_log10() +
    ylab("log(MSE)") +
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
          axis.title.y = element_text(size = 7))

}

# creating Standardization function
standardize = function(x){
  z <- (x - mean(x)) / sd(x)
  return( z)
}

# Plot digits in MNIST

plot_digit_img <- function(digit_df, palette, title_text, hjust, ncol = 4) {

  ggplot(data = digit_df, aes(x, y, fill = value)) +
    geom_tile() +
    facet_wrap(~ instance, ncol = ncol) +
    ggtitle(title_text) +
    coord_fixed() +
    scale_fill_continuous_sequential(palette = palette) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none",
          plot.title.position = "plot",
          plot.title = element_text(hjust = hjust,
                                    size = 15))

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


