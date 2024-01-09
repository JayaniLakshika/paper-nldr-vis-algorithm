calculate_effective_x_bins <- function(.data, x = UMAP1, cell_area = 1){

  if (any(is.na(.data$x))) {
    stop("NAs present")
  }

  if (any(is.infinite(.data$x))) {
    stop("Inf present")
  }

  if ((cell_area <= 0) || (is.infinite(cell_area))) {
    stop("Invalid cell area value")

  }

  cell_diameter <- sqrt(2 * cell_area / sqrt(3))

  xwidth <- diff(range(.data |>
                         dplyr::pull({{ x }})))

  num_bins <- ceiling(xwidth/cell_diameter)
  num_bins

}

calculate_effective_shape_value <- function(.data, x = UMAP1, y = UMAP2){

  if (any(is.na(.data$x)) || any(is.na(.data$y))) {
    stop("NAs present")
  }

  if (any(is.infinite(.data$x)) || any(is.infinite(.data$y))) {
    stop("Inf present")
  }

  if ((length(.data$x) == 1) || (length(.data$y) == 1)) {
    stop("Presence one observation only")

  }

  xwidth <- diff(range(.data |> dplyr::pull({{ x }})))
  yheight <- diff(range(.data |> dplyr::pull({{ y }})))


  shape <- yheight/xwidth  # Here, yheight is the range of y and xwidth is the range of x
  shape
}

extract_hexbin_centroids <- function(nldr_df, num_bins, shape_val = 1, x = UMAP1, y = UMAP2) {

  hb_data <- hexbin::hexbin(x = nldr_df |> dplyr::pull({{ x }}),
                            y = nldr_df |> dplyr::pull({{ y }}),
                            xbins = num_bins, IDs = TRUE,
                            shape = shape_val)

  hexdf_data <- tibble::tibble(tibble::as_tibble(hexbin::hcell2xy(hb_data)),  hexID = hb_data@cell, counts = hb_data@count, std_counts = hb_data@count/max(hb_data@count))

  return(list(hexdf_data = hexdf_data, hb_data = hb_data))
}

generate_full_grid_centroids <- function(hexdf_data){

  ## Generate initial grid
  full_centroids1 <- tibble::as_tibble(expand.grid(x = seq(min(hexdf_data$x),max(hexdf_data$x), ggplot2::resolution(hexdf_data$x, FALSE) * 2), y = seq(min(hexdf_data$y),max(hexdf_data$y), ggplot2::resolution(hexdf_data$y, FALSE) * 2)))

  ## Generate shifted grid
  full_centroids2 <- tibble::tibble(x = full_centroids1$x + ggplot2::resolution(hexdf_data$x, FALSE), y = full_centroids1$y + ggplot2::resolution(hexdf_data$y, FALSE))
  full_centroids <- dplyr::bind_rows(full_centroids1, full_centroids2)

  return(full_centroids)


}

full_hex_grid <- function(hexdf_data){

  dx <- ggplot2::resolution(hexdf_data$x, FALSE)
  dy <- ggplot2::resolution(hexdf_data$y, FALSE) / sqrt(3) / 2 * 1.15

  hexC <- hexbin::hexcoords(dx, dy, n = 1)

  n <- length(hexdf_data$x)

  size <- rep(1, length(hexdf_data$x))

  full_hex_coords <- tibble::tibble( x = rep.int(hexC$x, n) * rep(size, each = 6) + rep(hexdf_data$x, each = 6),
                                     y = rep.int(hexC$y, n) * rep(size, each = 6) + rep(hexdf_data$y, each = 6), id = rep(1:length(hexdf_data$x), each = 6))

  return(full_hex_coords)


}

map_polygon_id <- function(full_grid_with_hexbin_id, hex_grid){

  #browser()

  full_grid_with_polygon_id <- data.frame(matrix(ncol = 0, nrow = 0))

  for (i in 1:length(unique(full_grid_with_hexbin_id$hexID))) {

    full_grid_with_hexbin_id_filtered <- full_grid_with_hexbin_id |>
      dplyr::filter(hexID == unique(full_grid_with_hexbin_id$hexID)[i])

    for (j in 1:length(unique(hex_grid$id))) {

      hex_grid_filtered <- hex_grid |>
        dplyr::filter(id == unique(hex_grid$id)[j])

      status_in_x_range <- between(full_grid_with_hexbin_id_filtered$c_x, min(hex_grid_filtered$x), max(hex_grid_filtered$x))
      status_in_y_range <- between(full_grid_with_hexbin_id_filtered$c_y, min(hex_grid_filtered$y), max(hex_grid_filtered$y))

      if ((isTRUE(status_in_x_range)) & (isTRUE(status_in_y_range))) {

        full_grid_with_hexbin_id_filtered <- full_grid_with_hexbin_id_filtered |>
          dplyr::mutate(polygon_id = j)

        full_grid_with_polygon_id <- dplyr::bind_rows(full_grid_with_polygon_id, full_grid_with_hexbin_id_filtered)

      } else {

      }

    }

  }

  return(full_grid_with_polygon_id)

}

avg_highD_data <- function(.data, column_start_text = "x") {
  df_b <- .data |>
    dplyr::select(starts_with(column_start_text), hb_id) |>
    dplyr::group_by(hb_id) |>
    dplyr::summarise(across(everything(), mean))

  return(df_b)
}

triangulate_bin_centroids <- function(.data, x, y){
  tr1 <- tripack::tri.mesh(.data |> dplyr::pull({{ x }}), .data |> dplyr::pull({{ y }}))
  return(tr1)
}

generate_edge_info <- function(triangular_object) {
  tr_df <- tibble::tibble(x = triangular_object$x, y = triangular_object$y)  ## Create a dataframe with tri.mesh x and y coordinate values
  tr_df <- tr_df |>
    dplyr::mutate(ID = dplyr::row_number())  ## To add ID numbers, beacuse to join with from and to points in tri$arcs

  trang <- tripack::triangles(triangular_object)
  trang <- tibble::as_tibble(trang)

  tr_arcs_df1 <- tibble::tibble(from = trang$node1, to = trang$node2)  ## Create dataframe with from and to edges
  tr_arcs_df2 <- tibble::tibble(from = trang$node1, to = trang$node3)
  tr_arcs_df3 <- tibble::tibble(from = trang$node2, to = trang$node3)
  tr_arcs_df <- dplyr::bind_rows(tr_arcs_df1, tr_arcs_df2, tr_arcs_df3)  ## Create dataframe with from and to edges

  ## To obtain x and values of from to in a dataframe
  vec <- stats::setNames(rep("", 6), c("from", "to", "x_from", "y_from",
                                       "x_to", "y_to"))  ## Define column names
  # Initialize an empty dataframe to store data in a specific
  # format
  tr_from_to_df_coord <- dplyr::bind_rows(vec)[0, ]
  tr_from_to_df_coord <- tr_from_to_df_coord |>
    dplyr::mutate_if(is.character, as.numeric)

  for (i in 1:NROW(tr_arcs_df)) {
    from_row <- tr_df |>
      dplyr::filter(dplyr::row_number() == (tr_arcs_df |>
                                              dplyr::pull(from) |>
                                              dplyr::nth(i)))
    to_row <- tr_df |>
      dplyr::filter(dplyr::row_number() == (tr_arcs_df |>
                                              dplyr::pull(to) |>
                                              dplyr::nth(i)))
    tr_from_to_df_coord <- tr_from_to_df_coord |>
      tibble::add_row(from = from_row |>
                        dplyr::pull(ID), to = to_row |>
                        dplyr::pull(ID), x_from = from_row |>
                        dplyr::pull(x), y_from = from_row |>
                        dplyr::pull(y), x_to = to_row |>
                        dplyr::pull(x), y_to = to_row |>
                        dplyr::pull(y))  ## Add vector as an       appending row to the dataframe
  }

  return(tr_from_to_df_coord)
}

# cal_2D_dist <- function(.data){
#
#   .data$distance <- lapply(seq(nrow(.data)), function(x) {
#     start <- unlist(.data[x, c("x_from","y_from")])
#     end <- unlist(.data[x, c("x_to","y_to")])
#     sqrt(sum((start - end)^2))})
#
#   distance_df <- .data %>%
#     dplyr::select("from", "to", "distance")
#
#   distance_df$distance <- unlist(distance_df$distance)
#   return(distance_df)
# }


stat_trimesh <- function(mapping = NULL, data = NULL, geom = GeomTrimesh$default_aes(),
                         position = "identity", show.legend = NA, outliers = TRUE, inherit.aes = TRUE,
                         ...) {
  ggplot2::layer(
    stat = StatTrimesh,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(outliers = outliers, ...)
  )
}

StatTrimesh <- ggplot2::ggproto(
  "StatTrimesh",
  ggplot2::Stat,
  compute_group = function(data, scales, outliers = TRUE) {
    tr1 <- tripack::tri.mesh(data$x, data$y, duplicate = "remove")
    tr_df <- tibble::tibble(x = tr1$x, y = tr1$y)  ## Create a dataframe with tri.mesh x and y coordinate values
    tr_df <- tr_df |>
      dplyr::mutate(ID = dplyr::row_number())  ## To add ID numbers, because to join with from and to points in tri$arcs

    trang <- tripack::triangles(tr1)
    trang <- tibble::as_tibble(trang)

    tr_arcs_df1 <- tibble::tibble(from = trang$node1, to = trang$node2)  ## Create dataframe with from and to edges
    tr_arcs_df2 <- tibble::tibble(from = trang$node1, to = trang$node3)
    tr_arcs_df3 <- tibble::tibble(from = trang$node2, to = trang$node3)
    tr_arcs_df <- dplyr::bind_rows(tr_arcs_df1, tr_arcs_df2, tr_arcs_df3)  ## Create dataframe with from and to edges

    ## To obtain x and values of from to in a dataframe
    vec <- stats::setNames(rep("", 6), c("from", "to", "x_from", "y_from",
                                         "x_to", "y_to"))  ## Define column names
    # Initialize an empty dataframe to store data in a specific
    # format
    tr_from_to_df_coord <- dplyr::bind_rows(vec)[0, ]
    tr_from_to_df_coord <- tr_from_to_df_coord |>
      dplyr::mutate_if(is.character, as.numeric)

    for (i in 1:NROW(tr_arcs_df)) {
      from_row <- tr_df |>
        dplyr::filter(dplyr::row_number() == (tr_arcs_df |>
                                                dplyr::pull(from) |>
                                                dplyr::nth(i)))
      to_row <- tr_df |>
        dplyr::filter(dplyr::row_number() == (tr_arcs_df |>
                                                dplyr::pull(to) |>
                                                dplyr::nth(i)))
      tr_from_to_df_coord <- tr_from_to_df_coord |>
        tibble::add_row(from = from_row |>
                          dplyr::pull(ID),
                        to = to_row |>
                          dplyr::pull(ID),
                        x_from = from_row |>
                          dplyr::pull(x),
                        y_from = from_row |>
                          dplyr::pull(y),
                        x_to = to_row |>
                          dplyr::pull(x),
                        y_to = to_row |>
                          dplyr::pull(y))  ## Add vector as an appending row to the dataframe
    }

    trimesh <- tibble::tibble(x = tr_from_to_df_coord$x_from,
                              y = tr_from_to_df_coord$y_from,
                              xend = tr_from_to_df_coord$x_to,
                              yend = tr_from_to_df_coord$y_to,
                              PANEL = as.factor(rep(1, nrow(tr_from_to_df_coord))),
                              group = rep(-1, nrow(tr_from_to_df_coord)),
                              size = rep(0.5, nrow(tr_from_to_df_coord)),
                              linetype = rep(1, nrow(tr_from_to_df_coord)),
                              linewidth = rep(0.5, nrow(tr_from_to_df_coord)),
                              alpha = rep(NA, nrow(tr_from_to_df_coord)),
                              colour = rep("black", nrow(tr_from_to_df_coord)))
    trimesh
  },
  required_aes = c("x", "y")
)


geom_trimesh <- function(mapping = NULL, data = NULL, stat = "trimesh",
                         position = "identity", show.legend = NA, na.rm = FALSE, inherit.aes = TRUE,
                         ...) {
  ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = GeomTrimesh,
                 position = position, show.legend = show.legend, inherit.aes = inherit.aes,
                 params = list(na.rm = na.rm, ...))
}

GeomTrimesh <- ggplot2::ggproto("GeomTrimesh", ggplot2::Geom, required_aes = c("x", "y", "xend", "yend"),
                                default_aes = ggplot2::aes(shape = 19, linetype = 1, linewidth = 0.5,
                                                           size = 0.5, alpha = NA, colour = "black"),
                                draw_key = ggplot2::draw_key_point,
                                draw_panel = function(data, panel_scales, coord) {

                                  vertices <- tibble::tibble(x = data$x, y = data$y, colour = data$colour,
                                                             shape = data$shape, size = rep(2, nrow(data)), fill = rep("black",
                                                                                                                       nrow(data)), alpha = data$alpha, stroke = 0.5, stringsAsFactors = FALSE)

                                  trimesh <- tibble::tibble(x = data$x, xend = data$xend, y = data$y,
                                                            yend = data$yend, PANEL = data$PANEL, group = data$group, size = data$size,
                                                            linetype = data$linetype, linewidth = data$linewidth, alpha = data$alpha,
                                                            colour = data$colour)

                                  ggplot2:::ggname("geom_trimesh", grid::grobTree(ggplot2::GeomPoint$draw_panel(vertices,
                                                                                                                panel_scales, coord), ggplot2::GeomSegment$draw_panel(trimesh,
                                                                                                                                                                      panel_scales, coord)))
                                })



find_benchmark_value <- function(.data, distance_col) {
  #browser()

  .data <- .data |>
    dplyr::mutate(dplyr::across({
      {
        distance_col
      }
    }, \(x) round(x, 1)))


  sorted_distance_df <- .data |>
    dplyr::arrange({
      {
        distance_col
      }
    })  ## Sort the distances

  # b <- sorted_distance_df %>%
  #   group_by(distance) %>%
  #   summarise(n = n())
  #
  # benchmark_value <- b$distance[which(b$n == median(b$n))[1]]

  unique_dist <- sorted_distance_df |>
    dplyr::pull({
      {
        distance_col
      }
    }) |>
    unique()  ## Get the unique distances

  dist_u <- tibble::tibble(unique_dist = unique_dist)
  dist_u <- dplyr::bind_cols(dist_u, rbind(NA, apply(dist_u, 2, diff)), .name_repair = "unique_quiet")  ## Calculate differences between unique distance
  names(dist_u)[2] <- "difference"

  dist_u <- dist_u |>
    dplyr::mutate(dplyr::across(difference, \(x) round(x, 4)))  ## For simplicity

  dist_u[is.na(dist_u)] <- 0  ## To replace missing values with zero

  benchmark_value_vec <- c()

  ## To find the first largest difference (Define a benchmark value
  ## to remove long edges)
  for (i in 1:dim(dist_u)[1]) {
    if(!is.na(dist_u$difference[i + 1])){
      if (dist_u$difference[i] > dist_u$difference[i + 1]) {
        if (!(is.na(dist_u$difference[i]))) {
          benchmark_value_vec[i] <- dist_u$difference[i]
          break
        }
      }
    }
  }

  benchmark_value_df <- dist_u[which(dist_u$difference == benchmark_value_vec[!(is.na(benchmark_value_vec))]),
                               1]  # To get the first value which contain large difference
  names(benchmark_value_df) <- "unique_dist"
  benchmark_value <- benchmark_value_df |>
    dplyr::pull(unique_dist) |>
    dplyr::nth(1)
  benchmark_value

}


colour_long_edges <- function(.data, benchmark_value, triangular_object,
                              distance_col) {
  tr_df <- tibble::tibble(x = triangular_object$x, y = triangular_object$y)

  tr_from_to_df_coord <- generate_edge_info(triangular_object)

  distance_df_small_edges <- .data |>
    dplyr::filter({
      {
        distance_col
      }
    } < benchmark_value)
  distance_df_long_edges <- .data |>
    dplyr::filter({
      {
        distance_col
      }
    } >= benchmark_value)

  distance_df_small_edges <- distance_df_small_edges |>
    dplyr::mutate(type = "small_edges")

  distance_df_long_edges <- distance_df_long_edges |>
    dplyr::mutate(type = "long_edges")

  distance_edges <- dplyr::bind_rows(distance_df_small_edges, distance_df_long_edges)

  tr_from_to_df_coord_with_group <- merge(tr_from_to_df_coord, distance_edges,
                                          by = c("from", "to"))


  ## To draw the tri.mesh plot using ggplot
  tri_mesh_plot <- ggplot2::ggplot(tr_df, aes(x = x, y = y)) + ggplot2::geom_segment(aes(x = x_from,
                                                                                         y = y_from, xend = x_to, yend = y_to, color = type), data = tr_from_to_df_coord_with_group) +
    ggplot2::geom_point(size = 1) + ggplot2::coord_equal() + ggplot2::scale_colour_manual(values = c("#de2d26",
                                                                                                     "#636363"))  + ggplot2::labs(color=NULL) + ggplot2::theme(legend.position="bottom")
  return(tri_mesh_plot)

}

remove_long_edges <- function(.data, benchmark_value, triangular_object,
                              distance_col) {
  tr_df <- tibble::tibble(x = triangular_object$x, y = triangular_object$y)

  tr_from_to_df_coord <- generate_edge_info(triangular_object)

  distance_df_small_edges <- .data |>
    dplyr::filter({
      {
        distance_col
      }
    } < benchmark_value)

  tr_from_to_df_coord_with_group <- merge(tr_from_to_df_coord, distance_df_small_edges,
                                          by = c("from", "to"))


  ## To draw the tri.mesh plot using ggplot
  tri_mesh_plot <- ggplot2::ggplot(tr_df, aes(x = x, y = y)) + ggplot2::geom_segment(aes(x = x_from,
                                                                                         y = y_from, xend = x_to, yend = y_to), data = tr_from_to_df_coord_with_group) +
    ggplot2::geom_point(size = 1) + ggplot2::coord_equal() + ggplot2::labs(color=NULL)
  return(tri_mesh_plot)

}


show_langevitour <- function(df, df_b, df_b_with_center_data, benchmark_value = NA, distance_df, distance_col, min_points_threshold = NA, col_start = "x"){

  ### Define type column
  df <- df |>
    dplyr::select(tidyselect::starts_with(col_start)) |>
    dplyr::mutate(type = "data") ## original dataset

  df_b <- df_b |>
    dplyr::filter(hb_id %in% df_b_with_center_data$hexID) |>
    dplyr::select(-hb_id) |>
    dplyr::mutate(type = "model") ## Data with summarized mean

  df_exe <- dplyr::bind_rows(df_b, df)


  if((is.na(benchmark_value)) && (is.na(min_points_threshold))){

    tr1 <- triangulate_bin_centroids(df_b_with_center_data, x, y)
    tr_from_to_df <- generate_edge_info(triangular_object = tr1)

    langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = tr_from_to_df$from , lineTo = tr_from_to_df$to, group = df_exe$type)
  } else if ((!(is.na(benchmark_value))) && (is.na(min_points_threshold))) {
    ## Set the maximum difference as the criteria
    distance_df_small_edges <- distance_df %>%
      dplyr::filter({{ distance_col }} < benchmark_value)
    ## Since erase brushing is considerd.

    langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from, lineTo = distance_df_small_edges$to, group = df_exe$type)

  } else if ((is.na(benchmark_value)) && (!(is.na(min_points_threshold)))) {
    df_bin_centroids_filterd <- df_bin_centroids %>%
      dplyr::filter(counts > min_points_threshold)

    tr1 <- triangulate_bin_centroids(df_bin_centroids_filterd, x, y)
    tr_from_to_df <- generate_edge_info(triangular_object = tr1)

    langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = tr_from_to_df$from , lineTo = tr_from_to_df$to, group = df_exe$type)

  }  else if ((!(is.na(benchmark_value))) && (!(is.na(min_points_threshold)))) {

    df_bin_centroids_filterd <- df_bin_centroids %>%
      dplyr::filter(Cell_count > min_points_threshold)

    tr1 <- triangulate_bin_centroids(df_bin_centroids_filterd)
    tr_from_to_df <- generate_edge_info(triangular_object = tr1)

    distance_d <- cal_2D_dist(.data = tr_from_to_df)
    ## Set the maximum difference as the criteria
    distance_df_small_edges <- distance_d %>%
      dplyr::filter(distance < benchmark_value)
    ## Since erase brushing is considerd.

    langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from, lineTo = distance_df_small_edges$to, group = df_exe$type)

  } else {

  }


}

predict_hex_id <- function(training_data, nldr_df, nldr_df_test, num_bins, shape_val) {

  ## To extract bin centroids
  hexbin_data_object <-extract_hexbin_centroids(nldr_df, num_bins, shape_val)

  df_bin_centroids <- hexbin_data_object$hexdf_data

  UMAP_data_with_hb_id <- nldr_df |>
    dplyr::mutate(hb_id = hexbin_data_object$hb_data@cID)

  ## To generate a data set with high-D and 2D training data
  df_all <- dplyr::bind_cols(training_data |> dplyr::select(-ID), UMAP_data_with_hb_id)

  ## Averaged on high-D
  df_bin <- avg_highD_data(.data = df_all)

  train_hb_df <- df_bin_centroids |>
    dplyr::select(x, y, hexID)

  pred_hb_id <- class::knn(train_hb_df |> dplyr::select(-hexID), nldr_df_test |> dplyr::select(UMAP1, UMAP2), cl = train_hb_df$hexID)

  pred_data <- nldr_df_test |>
    dplyr::mutate(pred_hb_id = as.numeric(as.character(pred_hb_id)))

  return(list(pred_data = pred_data, df_bin_centroids = df_bin_centroids, df_bin = df_bin))

}

compute_aic <- function(p, total, num_bins, num_obs) {
  mse <- mean(total) / p
  aic <- 2*num_bins*p + num_obs*p*log(mse)
  return(aic)
}

generate_eval_df <- function(data, prediction_df, df_bin_centroids, df_bin, num_bins) {


  ## Generate all possible bin centroids in the full grid
  full_centroid_df <- generate_full_grid_centroids(df_bin_centroids)

  df_bin_centroids_filtered <- df_bin_centroids |>
    dplyr::select(hexID, x, y)

  ## To map centroid coordinates to predicted hexID
  prediction_df <- dplyr::inner_join(prediction_df, df_bin_centroids_filtered, by = c("pred_hb_id" = "hexID"))

  df_bin_train <- df_bin
  names(df_bin_train)[-1] <- paste0("avg_", names(df_bin_train)[-1])

  prediction_df <- prediction_df |>
    dplyr::left_join(df_bin_train, by = c("pred_hb_id" = "hb_id")) ## Map high-D averaged/weighted mean coordinates

  prediction_df <- prediction_df |>
    dplyr::left_join(data, by = c("ID" = "ID")) ## Map high-D data

  for (i in 1:(NCOL(df_bin_train) - 1)) {

    prediction_df[ , paste0("error_square_x", i)] <- (prediction_df[ , paste0("x", i)] - prediction_df[ , paste0("avg_x", i)])^2

  }

  prediction_df <- prediction_df |>
    dplyr::mutate(total = rowSums(dplyr::pick(tidyselect::starts_with("error_square_x"))))

  # prediction_df <- prediction_df |>
  #   dplyr::mutate(
  #     aic = compute_aic((NCOL(df_bin) - 1), prediction_df$total, NROW(full_centroid_df), NROW(prediction_df)),
  #     method2 = prediction_df$total * NROW(full_centroid_df)/NROW(prediction_df),
  #     method3 = prediction_df$total /NROW(full_centroid_df),
  #     mse = mean(prediction_df$total)/(NCOL(df_bin) - 1)
  #
  #   )

  # total_error <- sum(prediction_df$aic)
  # totol_error_method_2 <- sum(prediction_df$method2)
  # totol_error_method_3 <- sum(prediction_df$method3)
  # total_mse <- sum(prediction_df$mse)

  #number_of_bins: Total number of bins with empty bins
  eval_df <- tibble::tibble(number_of_bins = NROW(full_centroid_df), number_of_observations = NROW(prediction_df), total_error = compute_aic((NCOL(df_bin) - 1), prediction_df$total, NROW(full_centroid_df), NROW(prediction_df)), totol_error_method_2 = prediction_df$total * NROW(full_centroid_df)/NROW(prediction_df), totol_error_method_3 = prediction_df$total /NROW(full_centroid_df), total_mse = mean(prediction_df$total)/(NCOL(df_bin) - 1))

  return(eval_df)

}

generate_full_grid_info <- function(df_bin_centroids) {

  full_centroid_df <- generate_full_grid_centroids(df_bin_centroids)

  full_grid_with_hexbin_id <- map_hexbin_id(full_centroid_df, df_bin_centroids)

  ## Generate all coordinates of hexagons
  hex_grid <- full_hex_grid(full_centroid_df)

  full_grid_with_polygon_id_df <- map_polygon_id(full_grid_with_hexbin_id, hex_grid)

  full_grid_with_hexbin_id_rep <- full_grid_with_polygon_id_df |>
    dplyr::slice(rep(1:n(), each = 6)) |>
    dplyr::arrange(polygon_id)

  hex_full_count_df <- dplyr::bind_cols(hex_grid, full_grid_with_hexbin_id_rep)

  return(hex_full_count_df)

}

map_polygon_id <- function(full_grid_with_hexbin_id, hex_grid) {

  full_grid_with_polygon_id <- data.frame(matrix(ncol = 0, nrow = 0))

  for (i in 1:length(unique(full_grid_with_hexbin_id$hexID))) {

    full_grid_with_hexbin_id_filtered <- full_grid_with_hexbin_id |>
      dplyr::filter(hexID == unique(full_grid_with_hexbin_id$hexID)[i])

    for (j in 1:length(unique(hex_grid$id))) {

      hex_grid_filtered <- hex_grid |>
        dplyr::filter(id == unique(hex_grid$id)[j])

      status_in_x_range <- dplyr::between(full_grid_with_hexbin_id_filtered$c_x, min(hex_grid_filtered$x), max(hex_grid_filtered$x))
      status_in_y_range <- dplyr::between(full_grid_with_hexbin_id_filtered$c_y, min(hex_grid_filtered$y), max(hex_grid_filtered$y))

      if (any(status_in_x_range) & any(status_in_y_range)) {

        full_grid_with_hexbin_id_filtered <- full_grid_with_hexbin_id_filtered |>
          dplyr::mutate(polygon_id = j)

        full_grid_with_polygon_id <- dplyr::bind_rows(full_grid_with_polygon_id, full_grid_with_hexbin_id_filtered)
      }
    }
  }

  return(full_grid_with_polygon_id)
}

map_hexbin_id <- function(full_centroid_df, df_bin_centroids) {

  vec1 <- stats::setNames(rep("", 2), c("x", "y"))  ## Define column names

  full_grid_with_hexbin_id <- dplyr::bind_rows(vec1)[0, ]
  full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
    dplyr::mutate_if(is.character, as.numeric)

  for(i in 1:length(sort(unique(full_centroid_df$y)))) {

    ## Filter the data set with a specific y value
    specific_y_val_df <- full_centroid_df |>
      dplyr::filter(y == sort(unique(full_centroid_df$y))[i])

    ordered_x_df <- specific_y_val_df |>
      dplyr::arrange(x)

    full_grid_with_hexbin_id <- dplyr::bind_rows(full_grid_with_hexbin_id, ordered_x_df)
  }

  full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
    dplyr::mutate(hexID = row_number())

  full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
    dplyr::rename("c_x" = "x",
                  "c_y" = "y")

  full_grid_with_hexbin_id <- dplyr::full_join(full_grid_with_hexbin_id, df_bin_centroids, by = c("hexID" = "hexID")) |>
    dplyr::select(-c(x, y))

  full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
    dplyr::mutate(std_counts = counts/max(counts, na.rm = TRUE))

  return(full_grid_with_hexbin_id)
}

generate_full_grid_centroids <- function(hexdf_data){

  ## Generate initial grid
  full_centroids1 <- tibble::as_tibble(expand.grid(x = seq(min(hexdf_data$x),max(hexdf_data$x), ggplot2::resolution(hexdf_data$x, FALSE) * 2), y = seq(min(hexdf_data$y),max(hexdf_data$y), ggplot2::resolution(hexdf_data$y, FALSE) * 2)))

  ## Generate shifted grid
  full_centroids2 <- tibble::tibble(x = full_centroids1$x + ggplot2::resolution(hexdf_data$x, FALSE), y = full_centroids1$y + ggplot2::resolution(hexdf_data$y, FALSE))
  full_centroids <- dplyr::bind_rows(full_centroids1, full_centroids2)

  return(full_centroids)


}

full_hex_grid <- function(hexdf_data){

  dx <- ggplot2::resolution(hexdf_data$x, FALSE)
  dy <- ggplot2::resolution(hexdf_data$y, FALSE) / sqrt(3) / 2 * 1.15

  hexC <- hexbin::hexcoords(dx, dy, n = 1)

  n <- length(hexdf_data$x)

  size <- rep(1, length(hexdf_data$x))

  full_hex_coords <- tibble::tibble( x = rep.int(hexC$x, n) * rep(size, each = 6) + rep(hexdf_data$x, each = 6),
                                     y = rep.int(hexC$y, n) * rep(size, each = 6) + rep(hexdf_data$y, each = 6), id = rep(1:length(hexdf_data$x), each = 6))

  return(full_hex_coords)


}

find_low_density_hexagons <- function(df_bin_centroids, num_bins_x, benchmark_rm_hex = NA) {

  df_bin_centroids <- df_bin_centroids |>
    dplyr::mutate(ID = row_number())

  # To store mean densities of hexagons
  mean_density_vec <- c()

  for (i in 1:length(df_bin_centroids$hexID)) {

    df_bin_centroids_coordinates_spec_bin <- df_bin_centroids |>
      filter(hexID == df_bin_centroids$hexID[i])

    available_near_check <- df_bin_centroids |>
      dplyr::filter((hexID == (df_bin_centroids$hexID[i] + 1)) | (hexID == (df_bin_centroids$hexID[i] - 1))) |>
      head(1)

    if (NROW(available_near_check) == 0) {

      df_bin_centroids_coordinates_spec_bin_near1 <- df_bin_centroids |>
        filter((hexID == (df_bin_centroids$hexID[i] + (num_bins_x + 1))) | (hexID == (df_bin_centroids$hexID[i] + num_bins_x)) | (hexID == (df_bin_centroids$hexID[i] - (num_bins_x + 1))) | (hexID == (df_bin_centroids$hexID[i] - num_bins_x))) |>
        head(1)

    } else {

      df_bin_centroids_coordinates_spec_bin_near1 <- df_bin_centroids |>
        filter((hexID == (df_bin_centroids$hexID[i] + 1)) | (hexID == (df_bin_centroids$hexID[i] - 1))) |>
        head(1)

    }

    near_df_1 <- dplyr::bind_rows(df_bin_centroids_coordinates_spec_bin, df_bin_centroids_coordinates_spec_bin_near1)

    start <- unlist(near_df_1[1, c("x","y")])
    end <- unlist(near_df_1[2, c("x","y")])
    nearest_dist <- sqrt(sum((start - end)^2)) # Distance to nearest centroid

    df_bin_centroids$distance <- lapply(seq(nrow(df_bin_centroids)), function(x) {
      start <- unlist(df_bin_centroids[(df_bin_centroids_coordinates_spec_bin |> pull(ID)), c("x","y")])
      end <- unlist(df_bin_centroids[x, c("x","y")])
      sqrt(sum((start - end)^2))})

    df_bin_centroids <- df_bin_centroids %>%
      dplyr::select(names(df_bin_centroids), "distance")

    df_bin_centroids$distance <- round(unlist(df_bin_centroids$distance), 7)

    neighbor_df <- df_bin_centroids |>
      filter(distance == round(nearest_dist, 7))

    mean_density <- neighbor_df |>
      pull(std_counts) |>
      sum()/6 ## The reason to take the mean is to check the density in a considerable amount

    mean_density_vec <- append(mean_density_vec, mean_density)

  }

  df_bin_centroids <- df_bin_centroids |>
    dplyr::mutate(mean_density = mean_density_vec)

  remove_bins <- c()
  keep_bins <- c()

  for (i in 1:length(df_bin_centroids$hexID)) {

    df_bin_centroids_coordinates_spec_bin <- df_bin_centroids |>
      filter(hexID == df_bin_centroids$hexID[i])

    bin_ID <- df_bin_centroids_coordinates_spec_bin |>
      pull(hexID)

    if (is.na(benchmark_rm_hex)) {

      benchmark_rm_hex <- stats::quantile(mean_density_vec, probs = c(0,0.25,0.5,0.75,1))[2]

    }

    if(df_bin_centroids_coordinates_spec_bin$mean_density < benchmark_rm_hex){
      remove_bins <- append(remove_bins, bin_ID)
    } else {
      keep_bins <- append(keep_bins, bin_ID)
    }
  }

  return(remove_bins)
}


extract_hexbin_mean <- function(nldr_df, num_bins, shape_val = 1, x = UMAP1, y = UMAP2) {

  hb_data <- hexbin::hexbin(x = nldr_df |> dplyr::pull({{ x }}),
                            y = nldr_df |> dplyr::pull({{ y }}),
                            xbins = num_bins, IDs = TRUE,
                            shape = shape_val)

  df_cell_data <- nldr_df |>
    dplyr::select(-ID) |>
    dplyr::mutate(hexID = hb_data@cID) |>
    dplyr::group_by(hexID) |>
    dplyr::summarise(dplyr::across(tidyselect::everything(), mean))

  names(df_cell_data) <- c("hexID", "x", "y")


  hexdf_data <- tibble::tibble(df_cell_data, counts = hb_data@count, std_counts = hb_data@count/max(hb_data@count))

  return(list(hexdf_data = hexdf_data, hb_data = hb_data))
}


cal_2D_dist_umap <- function(.data){

  .data$distance <- lapply(seq(nrow(.data)), function(x) {
    start <- unlist(.data[x, c("avg_umap1","avg_umap2")])
    end <- unlist(.data[x, c("UMAP1","UMAP2")])
    sqrt(sum((start - end)^2))})

  distance_df <- .data %>%
    dplyr::select("hb_id", "avg_umap1","avg_umap2", "UMAP1","UMAP2", "distance")

  distance_df$distance <- unlist(distance_df$distance)
  return(distance_df)
}


cal_2D_dist <- function(.data, start_x = "x_from", start_y = "y_from", end_x = "x_to", end_y = "y_to", select_col_vec = c("from", "to", "distance")) {
  # Calculate the 2D distances
  .data$distance <- lapply(seq(nrow(.data)), function(x) {
    start <- unlist(.data[x, c(start_x, start_y)])
    end <- unlist(.data[x, c(end_x, end_y)])
    sqrt(sum((start - end)^2))
  })

  # Create a data frame with the from-to relationships and distances
  distance_df <- .data |> dplyr::select(all_of(select_col_vec))

  # Convert the distances to a vector and return the data frame
  distance_df$distance <- unlist(distance_df$distance)
  return(distance_df)
}


compute_weights <- function(nldr_df, hb_object) {

  ## To get the average of each bin
  bin_val_hexagons <- nldr_df |>
    dplyr::mutate(hb_id = hb_object@cID) |>
    dplyr::select(-ID) |>
    dplyr::group_by(hb_id) |>
    dplyr::summarise(dplyr::across(tidyselect::everything(), mean))

  new_col <- paste0("avg_", names(UMAP_data)[1:2] |> tolower())

  names(bin_val_hexagons) <- append("hb_id", new_col)

  ## To calculate distances from average point

  umap_with_avg_all <- dplyr::inner_join(bin_val_hexagons , nldr_df |>
                                           dplyr::mutate(hb_id = hb_object@cID) |>
                                           dplyr::select(-ID), by = c("hb_id" = "hb_id"))


  umap_with_avg_all_split <- umap_with_avg_all |>
    dplyr::group_by(hb_id) |>
    dplyr::group_split()

  col_names1 <- append(names(bin_val_hexagons), (names(nldr_df) - 1))
  col_names <- append(col_names1, "distance")

  vec <- stats::setNames(1:6, col_names)
  weight_df <- dplyr::bind_rows(vec)[0, ]

  for(i in 1:length(umap_with_avg_all_split)){

    weighted_mean_df <- umap_with_avg_all_split[[i]] |> ## These are the weights for weighted mean
      cal_2D_dist_umap()

    weight_df <- dplyr::bind_rows(weight_df, weighted_mean_df)

  }

  return(weight_df)

}


weighted_highD_data <- function(.data, weight_df) {

  weighted_mean_all <- dplyr::inner_join(.data, weight_df, by = c("hb_id" = "hb_id", "UMAP1" = "UMAP1", "UMAP2" = "UMAP2")) |>
    mutate(distance_trans =  1/ (distance + 0.05))

  weighted_mean_df_list <- list()

  for (j in 1:(NCOL(weighted_mean_all) - 8)) {

    weighted_mean_df_list[[j]] <- weighted_mean_all |>
      dplyr::select(hb_id, names(weighted_mean_all)[-((length(weighted_mean_all)-7):length(weighted_mean_all))][j], distance_trans) |>
      dplyr::group_by(hb_id) |>
      dplyr::summarise(dplyr::across(names(weighted_mean_all)[-((length(weighted_mean_all)-7):length(weighted_mean_all))][j], ~ weighted.mean(., distance_trans)))

  }

  weighted_mean <- weighted_mean_df_list %>%
    Reduce(function(dtf1,dtf2) dplyr::full_join(dtf1,dtf2,by="hb_id"), .)


  ## Column names starts with x
  weighted_mean <- weighted_mean |>
    dplyr::select(hb_id, tidyselect::starts_with("x"))

  return(weighted_mean)
}
