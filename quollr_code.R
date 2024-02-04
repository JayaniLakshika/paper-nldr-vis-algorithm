calculate_effective_x_bins <- function(.data, x = UMAP1, cell_area = 1){

  if (any(is.na(.data |> dplyr::pull({{ x }})))) {
    stop("NAs present")
  }

  if (any(is.infinite(.data |> dplyr::pull({{ x }})))) {
    stop("Inf present")
  }

  if ((cell_area <= 0) || (is.infinite(cell_area))) {
    stop("Invalid cell area value")

  }

  ## To compute the diameter of the hexagon
  cell_diameter <- sqrt(2 * cell_area / sqrt(3))

  ## To compute the range along x-axis
  xwidth <- diff(range(.data |>
                         dplyr::pull({{ x }})))

  num_bins <- ceiling(xwidth/cell_diameter)
  #num_bins <- floor(xwidth/cell_diameter + 1.5001)
  num_bins

}

calculate_effective_shape_value <- function(.data, x = UMAP1, y = UMAP2){

  if (any(is.na(.data |> dplyr::pull({{ x }}))) || any(is.na(.data |> dplyr::pull({{ y }})))) {
    stop("NAs present")
  }

  if (any(is.infinite(.data |> dplyr::pull({{ x }}))) || any(is.infinite(.data |> dplyr::pull({{ y }})))) {
    stop("Inf present")
  }

  if ((length(.data |> dplyr::pull({{ x }})) == 1) || (length(.data |> dplyr::pull({{ y }})) == 1)) {
    stop("Presence one observation only")

  }

  ## To compute the range along x-axis
  xwidth <- diff(range(.data |> dplyr::pull({{ x }})))
  ## To compute the range along y-axis
  yheight <- diff(range(.data |> dplyr::pull({{ y }})))


  shape <- yheight/xwidth
  shape
}

extract_hexbin_centroids <- function(nldr_df, num_bins, shape_val = 1, x = UMAP1, y = UMAP2) {

  ## To create the hexbin object
  hb_data <- hexbin::hexbin(x = nldr_df |> dplyr::pull({{ x }}),
                            y = nldr_df |> dplyr::pull({{ y }}),
                            xbins = num_bins, IDs = TRUE,
                            shape = shape_val)

  ## To create the hexbin centroid info dataset
  hexdf_data <- tibble::tibble(tibble::as_tibble(hexbin::hcell2xy(hb_data)),  hexID = hb_data@cell, counts = hb_data@count, std_counts = hb_data@count/max(hb_data@count))

  return(list(hexdf_data = hexdf_data, hb_data = hb_data))
}

extract_hexbin_mean <- function(nldr_df, num_bins, shape_val = 1, x = UMAP1, y = UMAP2) {

  ## To create the hexbin object
  hb_data <- hexbin::hexbin(x = nldr_df |> dplyr::pull({{ x }}),
                            y = nldr_df |> dplyr::pull({{ y }}),
                            xbins = num_bins, IDs = TRUE,
                            shape = shape_val)

  ## To compute hexagonal bin means
  df_cell_data <- nldr_df |>
    dplyr::select(-ID) |>
    dplyr::mutate(hexID = hb_data@cID) |>
    dplyr::group_by(hexID) |>
    dplyr::summarise(dplyr::across(tidyselect::everything(), mean))

  ## To rename the columns
  names(df_cell_data) <- c("hexID", "x", "y")

  ## To create the hexbin means info dataset
  hexdf_data <- tibble::tibble(df_cell_data, counts = hb_data@count, std_counts = hb_data@count/max(hb_data@count))

  return(list(hexdf_data = hexdf_data, hb_data = hb_data))
}

triangulate_bin_centroids <- function(.data, x, y){
  tr1 <- tripack::tri.mesh(.data |> dplyr::pull({{ x }}), .data |> dplyr::pull({{ y }}))
  return(tr1)
}

generate_edge_info <- function(triangular_object) {
  # Create a data frame with x and y coordinate values from the triangular object
  tr_df <- tibble::tibble(x = triangular_object$x, y = triangular_object$y)
  tr_df <- tr_df |> dplyr::mutate(ID = dplyr::row_number())  # Add ID numbers for joining with from and to points in tr_arcs

  # Extract the triangles from the triangular object
  trang <- tripack::triangles(triangular_object)
  trang <- tibble::as_tibble(trang)

  # Create data frames with from-to edges
  tr_arcs_df1 <- tibble::tibble(from = trang$node1, to = trang$node2)
  tr_arcs_df2 <- tibble::tibble(from = trang$node1, to = trang$node3)
  tr_arcs_df3 <- tibble::tibble(from = trang$node2, to = trang$node3)
  tr_arcs_df <- dplyr::bind_rows(tr_arcs_df1, tr_arcs_df2, tr_arcs_df3)

  tr_arcs_df <- tr_arcs_df |> ## To remove duplicates
    dplyr::distinct()

  ## To extract unique combinations
  tr_arcs_df <- tr_arcs_df |>
    dplyr::mutate(x = pmin(from, to), y = pmax(from, to)) |>
    dplyr::distinct(x, y) |>
    dplyr::rename(c("from" = "x", "to" = "y"))

  # Create an empty data frame to store the edge information
  vec <- stats::setNames(rep("", 6), c("from", "to", "x_from", "y_from", "x_to", "y_to"))  # Define column names
  tr_from_to_df_coord <- dplyr::bind_rows(vec)[0, ]
  tr_from_to_df_coord <- tr_from_to_df_coord |> dplyr::mutate_if(is.character, as.numeric)

  # Generate the edge information
  for (i in 1:NROW(tr_arcs_df)) {
    from_row <- tr_df |> dplyr::filter(dplyr::row_number() == (tr_arcs_df |> dplyr::pull(from) |> dplyr::nth(i)))
    to_row <- tr_df |> dplyr::filter(dplyr::row_number() == (tr_arcs_df |> dplyr::pull(to) |> dplyr::nth(i)))
    tr_from_to_df_coord <- tr_from_to_df_coord |> tibble::add_row(
      from = from_row |> dplyr::pull(ID),
      to = to_row |> dplyr::pull(ID),
      x_from = from_row |> dplyr::pull(x),
      y_from = from_row |> dplyr::pull(y),
      x_to = to_row |> dplyr::pull(x),
      y_to = to_row |> dplyr::pull(y)
    )
  }

  return(tr_from_to_df_coord)
}

cal_2d_dist <- function(.data, start_x = "x_from", start_y = "y_from", end_x = "x_to",
                        end_y = "y_to", select_col_vec = c("from", "to", "distance")) {
  # Calculate the 2D distances
  .data$distance <- lapply(seq(nrow(.data)), function(x) {
    start <- unlist(.data[x, c(start_x, start_y)])
    end <- unlist(.data[x, c(end_x, end_y)])
    sqrt(sum((start - end)^2))
  })

  # Create a data frame with the from-to relationships and distances
  distance_df <- .data |> dplyr::select(tidyselect::all_of(select_col_vec))

  # Convert the distances to a vector and return the data frame
  distance_df$distance <- unlist(distance_df$distance)
  return(distance_df)
}

colour_long_edges <- function(.data, benchmark_value, triangular_object, distance_col) {
  # Create the tibble with x and y coordinates
  tr_df <- tibble::tibble(x = triangular_object$x, y = triangular_object$y)

  # Generate edge information
  tr_from_to_df_coord <- generate_edge_info(triangular_object)

  # Filter and label small and long edges
  distance_df_small_edges <- .data |>
    dplyr::filter({{ distance_col }} < benchmark_value) |>
    dplyr::mutate(type = "small_edges")

  distance_df_long_edges <- .data |>
    dplyr::filter({{ distance_col }} >= benchmark_value) |>
    dplyr::mutate(type = "long_edges")

  # Combine small and long edges
  distance_edges <- dplyr::bind_rows(distance_df_small_edges, distance_df_long_edges)

  # Merge edge information with distance data
  tr_from_to_df_coord_with_group <- merge(tr_from_to_df_coord, distance_edges, by = c("from", "to"))

  # Create the triangular mesh plot with colored long edges
  tri_mesh_plot <- ggplot2::ggplot(tr_df, aes(x = x, y = y)) +
    ggplot2::geom_segment(
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to, color = type),
      data = tr_from_to_df_coord_with_group
    ) +
    ggplot2::geom_point(size = 1, colour = "#33a02c") +
    ggplot2::coord_equal() +
    ggplot2::scale_colour_manual(values = c("#de2d26", "#636363"))

  return(tri_mesh_plot)
}

remove_long_edges <- function(.data, benchmark_value, triangular_object,
                              distance_col) {
  # Create the tibble with x and y coordinates
  tr_df <- tibble::tibble(x = triangular_object$x, y = triangular_object$y)

  # Generate edge information
  tr_from_to_df_coord <- generate_edge_info(triangular_object)

  # Filter small edges
  distance_df_small_edges <- .data |>
    dplyr::filter({
      {
        distance_col
      }
    } < benchmark_value)

  # Merge edge information with distance data
  tr_from_to_df_coord_with_group <- merge(tr_from_to_df_coord, distance_df_small_edges,
                                          by = c("from", "to"))


  ## Create the triangular mesh plot after removing the long edges
  tri_mesh_plot <- ggplot2::ggplot(tr_df, aes(x = x, y = y)) + ggplot2::geom_segment(aes(x = x_from,
                                                                                         y = y_from, xend = x_to, yend = y_to), data = tr_from_to_df_coord_with_group) +
    ggplot2::geom_point(size = 1, colour = "#33a02c") + ggplot2::coord_equal() + ggplot2::labs(color=NULL)
  return(tri_mesh_plot)

}

generate_full_grid_centroids <- function(hexdf_data){

  ## Generate initial grid
  full_centroids1 <- tibble::as_tibble(expand.grid(x = seq(min(hexdf_data$x),max(hexdf_data$x), ggplot2::resolution(hexdf_data$x, FALSE) * 2), y = seq(min(hexdf_data$y),max(hexdf_data$y), ggplot2::resolution(hexdf_data$y, FALSE) * 2)))

  ## Generate shifted grid
  full_centroids2 <- tibble::tibble(x = full_centroids1$x + ggplot2::resolution(hexdf_data$x, FALSE), y = full_centroids1$y + ggplot2::resolution(hexdf_data$y, FALSE))

  ## Combine all
  full_centroids <- dplyr::bind_rows(full_centroids1, full_centroids2)

  return(full_centroids)


}

full_hex_grid <- function(hexdf_data){

  ## Compute horizontal width of the hexagon
  dx <- ggplot2::resolution(hexdf_data$x, FALSE)

  ## Compute vertical width of the hexagon
  # Adjust for difference in width and height of regular hexagon. 1.15 adjusts
  # for the effect of the overlapping range in y-direction on the resolution
  dy <- ggplot2::resolution(hexdf_data$y, FALSE) / sqrt(3) / 2 * 1.15

  ## Obtain hexagon polygon coordinates
  hexC <- hexbin::hexcoords(dx, dy, n = 1) ## n: Number of hexagons repeat

  ## Obtain the number of hexagons in the full grid
  n <- length(hexdf_data$x)

  ## Generate the size vector of the hexagons (since regular hexagons)
  size <- rep(1, length(hexdf_data$x))

  ## Generate the coordinates for the hexagons
  full_hex_coords <- tibble::tibble( x = rep.int(hexC$x, n) * rep(size, each = 6) + rep(hexdf_data$x, each = 6),
                                     y = rep.int(hexC$y, n) * rep(size, each = 6) + rep(hexdf_data$y, each = 6),
                                     id = rep(1:length(hexdf_data$x), each = 6))

  return(full_hex_coords)


}

map_hexbin_id <- function(full_centroid_df, df_bin_centroids) {

  vec1 <- stats::setNames(rep("", 2), c("x", "y"))  ## Define column names

  ## Define a dataset to store all the centroids with the respective coordinates
  full_grid_with_hexbin_id <- dplyr::bind_rows(vec1)[0, ]
  full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
    dplyr::mutate_if(is.character, as.numeric)

  for(i in 1:length(sort(unique(full_centroid_df$y)))){

    ## Filter the data set with specific y value
    specific_y_val_df <- full_centroid_df |>
      dplyr::filter(y == sort(unique(full_centroid_df$y))[i])

    ## orderd the x values
    ordered_x_df <- specific_y_val_df |>
      dplyr::arrange(x)

    full_grid_with_hexbin_id <- dplyr::bind_rows(full_grid_with_hexbin_id, ordered_x_df)

  }

  ## Add the column with hexagonal bin ID
  full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
    dplyr::mutate(hexID = dplyr::row_number())

  ## Rename columns
  full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
    dplyr::rename("c_x" = "x",
                  "c_y" = "y")

  ## Join with centroid data set to extarct the count column
  full_grid_with_hexbin_id <- dplyr::full_join(full_grid_with_hexbin_id, df_bin_centroids, by = c("hexID" = "hexID")) |>
    dplyr::select(-c(x, y))

  ## Compute the standardise count
  full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
    dplyr::mutate(std_counts = counts/max(counts, na.rm = TRUE))

  return(full_grid_with_hexbin_id)
}

map_polygon_id <- function(full_grid_with_hexbin_id, hex_grid) {

  ## Define a dataset to store polygon id
  full_grid_with_polygon_id <- data.frame(matrix(ncol = 0, nrow = 0))

  for (i in 1:length(unique(full_grid_with_hexbin_id$hexID))) {

    ## Filter specific hexagon
    full_grid_with_hexbin_id_filtered <- full_grid_with_hexbin_id |>
      dplyr::filter(hexID == unique(full_grid_with_hexbin_id$hexID)[i])

    for (j in 1:length(unique(hex_grid$id))) {

      ## Filter sepcific polygon
      hex_grid_filtered <- hex_grid |>
        dplyr::filter(id == unique(hex_grid$id)[j])

      ## Check the centroid exists within the polygon
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

generate_full_grid_info <- function(df_bin_centroids) {

  ## generate all the centroids
  full_centroid_df <- generate_full_grid_centroids(df_bin_centroids)

  ## Map the hexgoanl ID to full centroid dataset
  full_grid_with_hexbin_id <- map_hexbin_id(full_centroid_df, df_bin_centroids)

  ## Generate all coordinates of hexagons
  hex_grid <- full_hex_grid(full_centroid_df)

  ## Map the polygon ID to the hexagon coordinates
  full_grid_with_polygon_id_df <- map_polygon_id(full_grid_with_hexbin_id, hex_grid)

  full_grid_with_hexbin_id_rep <- full_grid_with_polygon_id_df |>
    dplyr::slice(rep(1:dplyr::n(), each = 6)) |>
    dplyr::arrange(polygon_id)

  ## Generate the dataset with polygon, and hexagon bin centroid coordinates
  hex_full_count_df <- dplyr::bind_cols(hex_grid, full_grid_with_hexbin_id_rep)

  return(hex_full_count_df)

}

find_pts_in_hexbins <- function(full_grid_with_hexbin_id, nldr_data_with_hb_id) {

  ## Dataframe to store points info
  pts_df <- data.frame(matrix(ncol = 0, nrow = 0))

  for (i in 1:length(nldr_data_with_hb_id$hb_id)) {

    ## Filter a hexagon and find the point within that hexagon
    pts_vec <- nldr_data_with_hb_id |>
      dplyr::filter(hb_id == nldr_data_with_hb_id$hb_id[i]) |>
      dplyr::pull(ID)

    ## Store the hexagon ID with the respective points
    hb_pts <- tibble::tibble(hexID = nldr_data_with_hb_id$hb_id[i], pts = list(pts_vec))

    pts_df <- dplyr::bind_rows(pts_df, hb_pts)

  }

  return(pts_df)

}

find_non_empty_bins <- function(nldr_df, x = "UMAP1", y = "UMAP2", shape_val, non_empty_bins) {

  num_bins_x <- 1
  ## To extract bin centroids
  hexbin_data_object <- extract_hexbin_centroids(nldr_df = nldr_df,
                                                 num_bins = num_bins_x,
                                                 shape_val = shape_val, x = x, y = y)
  df_bin_centroids <- hexbin_data_object$hexdf_data

  num_of_non_empty_bins <- df_bin_centroids$hexID |> length()

  while (num_of_non_empty_bins < non_empty_bins) {

    num_bins_x <- num_bins_x + 1

    ## To extract bin centroids
    hexbin_data_object <- extract_hexbin_centroids(nldr_df = nldr_df,
                                                   num_bins = num_bins_x,
                                                   shape_val = shape_val, x = y, y = y)

    df_bin_centroids <- hexbin_data_object$hexdf_data

    num_of_non_empty_bins <- df_bin_centroids$hexID |> length()

    if (num_of_non_empty_bins >= non_empty_bins) {
      return(num_bins_x)
      break
    } else {
      next
    }
  }
}

avg_highD_data <- function(.data, column_start_text = "x") {
  df_b <- .data |>
    dplyr::select(rsample::starts_with(column_start_text), hb_id) |>
    dplyr::group_by(hb_id) |>
    dplyr::summarise(dplyr::across(tidyselect::everything(), mean))

  return(df_b)
}

compute_weights <- function(nldr_df, hb_object) {

  hb_object <- hb_object

  ## To get the average of each bin
  bin_val_hexagons <- nldr_df |>
    dplyr::mutate(hb_id = hb_object@cID) |>
    dplyr::group_by(hb_id) |>
    dplyr::summarise(dplyr::across(tidyselect::everything(), mean))

  new_col <- paste0("avg_", names(nldr_df)[1:2] |> tolower())

  names(bin_val_hexagons) <- append("hb_id", new_col)

  ## To calculate distances from average point

  nldr_with_avg_all <- dplyr::inner_join(bin_val_hexagons , nldr_df |>
                                           dplyr::mutate(hb_id = hb_object@cID),
                                         by = c("hb_id" = "hb_id"))


  nldr_with_avg_all_split <- nldr_with_avg_all |>
    dplyr::group_by(hb_id) |>
    dplyr::group_split()

  col_names1 <- append(names(bin_val_hexagons), names(nldr_df))
  col_names <- append(col_names1, "distance")

  vec <- stats::setNames(1:6, col_names)
  weight_df <- dplyr::bind_rows(vec)[0, ]

  for(i in 1:length(nldr_with_avg_all_split)){

    weighted_mean_df <- nldr_with_avg_all_split[[i]] |> ## These are the weights for weighted mean
      cal_2d_dist(start_x = new_col[1], start_y = new_col[2], end_x = names(nldr_df)[1], end_y = names(nldr_df)[2], select_col_vec = col_names)

    weight_df <- dplyr::bind_rows(weight_df, weighted_mean_df)

  }

  return(weight_df)

}

weighted_highD_data <- function(training_data, nldr_df_with_id, hb_object, column_start_text = "x") {

  nldr_df_with_id <- nldr_df_with_id |> dplyr::mutate(hb_id = hb_object@cID)
  df_all <- dplyr::bind_cols(training_data |> dplyr::select(-ID), nldr_df_with_id)

  weight_df <- compute_weights(nldr_df_with_id |> dplyr::select(-c(ID, hb_id)), hb_object)

  joined_col_names <- names(nldr_df_with_id)[1:2]

  weighted_mean_all <- dplyr::inner_join(df_all, weight_df, by = c("hb_id" = "hb_id",
                                                                   stats::setNames(joined_col_names, joined_col_names))) |>
    mutate(distance_trans =  1/ (distance + 0.05))

  weighted_mean_df_list <- list()

  for (j in 1:(NCOL(training_data |> dplyr::select(-ID)))) {

    weighted_mean_df_list[[j]] <- weighted_mean_all |>
      dplyr::select(hb_id, names(training_data |> dplyr::select(-ID))[j], distance_trans) |>
      dplyr::group_by(hb_id) |>
      dplyr::summarise(dplyr::across(names(training_data |> dplyr::select(-ID))[j], ~ weighted.mean(., distance_trans)))

  }

  weighted_mean <- Reduce(function(dtf1,dtf2) dplyr::full_join(dtf1,dtf2,by="hb_id"),
                          weighted_mean_df_list)


  ## Column names start with x
  weighted_mean <- weighted_mean |>
    dplyr::select(hb_id, tidyselect::starts_with(column_start_text))

  return(weighted_mean)
}

show_langevitour <- function(df, df_b, df_b_with_center_data, benchmark_value = NA,
                             distance_df, distance_col, min_points_threshold = NA) {

  ### Define type column
  df <- df |>
    dplyr::select(tidyselect::starts_with("x")) |>
    dplyr::mutate(type = "data") ## original dataset

  df_b <- df_b |>
    dplyr::filter(hb_id %in% df_b_with_center_data$hexID) |>
    dplyr::select(-hb_id) |>
    dplyr::mutate(type = "model") ## Data with summarized mean

  df_exe <- dplyr::bind_rows(df_b, df)


  if((is.na(benchmark_value)) && (is.na(min_points_threshold))){

    tr1 <- triangulate_bin_centroids(df_b_with_center_data, x, y)
    tr_from_to_df <- generate_edge_info(triangular_object = tr1)

    langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = tr_from_to_df$from , lineTo = tr_from_to_df$to, group = df_exe$type, pointSize = 3, levelColors = c("#6a3d9a", "#33a02c"))
  } else if ((!(is.na(benchmark_value))) && (is.na(min_points_threshold))) {
    ## Set the maximum difference as the criteria
    distance_df_small_edges <- distance_df |>
      dplyr::filter({{ distance_col }} < benchmark_value)
    ## Since erase brushing is considerd.

    langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from, lineTo = distance_df_small_edges$to, group = df_exe$type, pointSize = 3, levelColors = c("#6a3d9a", "#33a02c"))

  } else if ((is.na(benchmark_value)) && (!(is.na(min_points_threshold)))) {
    df_bin_centroids_filterd <- df_bin_centroids |>
      dplyr::filter(counts > min_points_threshold)

    tr1 <- triangulate_bin_centroids(df_bin_centroids_filterd, x, y)
    tr_from_to_df <- generate_edge_info(triangular_object = tr1)

    langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = tr_from_to_df$from , lineTo = tr_from_to_df$to, group = df_exe$type, pointSize = 3, levelColors = c("#6a3d9a", "#33a02c"))

  }  else if ((!(is.na(benchmark_value))) && (!(is.na(min_points_threshold)))) {

    df_bin_centroids_filterd <- df_bin_centroids |>
      dplyr::filter(Cell_count > min_points_threshold)

    tr1 <- triangulate_bin_centroids(df_bin_centroids_filterd)
    tr_from_to_df <- generate_edge_info(triangular_object = tr1)

    distance_d <- cal_2d_dist(.data = tr_from_to_df)
    ## Set the maximum difference as the criteria
    distance_df_small_edges <- distance_d |>
      dplyr::filter(distance < benchmark_value)
    ## Since erase brushing is considerd.

    langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from, lineTo = distance_df_small_edges$to, group = df_exe$type, pointSize = 3, levelColors = c("#6a3d9a", "#33a02c"))

  } else {

  }


}

fit_high_d_model <- function(training_data, nldr_df_with_id, x = "UMAP1",
                             y = "UMAP1", cell_area = 1, num_bins_x = NA, shape_val = NA,
                             is_bin_centroid = TRUE,
                             is_rm_lwd_hex = FALSE,
                             benchmark_to_rm_lwd_hex = NA,
                             is_avg_high_d = TRUE, column_start_text = "x") {

  ## If number of bins along the x-axis is not given
  if (is.na(num_bins_x)) {
    ## compute the number of bins along the x-axis
    num_bins_x <- calculate_effective_x_bins(.data = nldr_df_with_id, x = x,
                                             cell_area = cell_area)


  }

  ## If shape parameter is not given
  if (is.na(shape_val)) {
    ## compute shape parameter
    shape_val <- calculate_effective_shape_value(.data = nldr_df_with_id, x = x, y = y)

  }

  ## Do you need to use bin centroids or bin means?
  if (isTRUE(is_bin_centroid)) {
    ## For bin centroids
    hexbin_data_object <- extract_hexbin_centroids(nldr_df = nldr_df_with_id,
                                                   num_bins = num_bins_x,
                                                   shape_val = shape_val,
                                                   x = x, y = y)

  } else {
    ## For bin means
    hexbin_data_object <- extract_hexbin_mean(nldr_df = nldr_df_with_id,
                                              num_bins = num_bins_x,
                                              shape_val = shape_val,
                                              x = x, y = y)

  }



  ## Do you need to remove low density hexagons?
  if (isTRUE(is_rm_lwd_hex)) {
    ## extract bin centroid/bin mean info
    df_bin_centroids <- hexbin_data_object$hexdf_data

    ## if the benchmark value to remove low density hexagons is not provided
    if (is.na(benchmark_to_rm_lwd_hex)) {
      ## first quartile used as the default
      benchmark_to_rm_lwd_hex <- stats::quantile(df_bin_centroids$std_counts,
                                                 probs = c(0,0.25,0.5,0.75,1))[2]
    }

    ## To identify low density hexagons
    df_bin_centroids_low <- df_bin_centroids |>
      dplyr::filter(std_counts <= benchmark_to_rm_lwd_hex)

    ## To identify low-density hexagons needed to remove by investigating neighbouring mean density
    identify_rm_bins <- find_low_density_hexagons(df_bin_centroids_all = df_bin_centroids,
                                                  num_bins_x = num_bins_x,
                                                  df_bin_centroids_low = df_bin_centroids_low)

    ## To remove low-density hexagons
    df_bin_centroids <- df_bin_centroids |>
      filter(!(hexID %in% identify_rm_bins))



  } else {
    ### Witout removing low density hexagons
    ## extract bin centroid/bin mean info
    df_bin_centroids <- hexbin_data_object$hexdf_data

  }


  nldr_df_with_hb_id <- nldr_df_with_id |>
    dplyr::mutate(hb_id = hexbin_data_object$hb_data@cID)

  ## To generate a data set with high-D and 2D training data
  df_all <- dplyr::bind_cols(training_data |> dplyr::select(-ID), nldr_df_with_hb_id)

  ## Do you need to use bin centroids or bin means?
  if (isTRUE(is_avg_high_d)) {

    ## averaged high-D data
    df_bin <- avg_highD_data(.data = df_all, column_start_text = column_start_text)


  } else {

    ## weighted averaged high-D data
    df_bin <- weighted_highD_data(training_data = training_data,
                                  nldr_df_with_id = nldr_df_with_id,
                                  hb_object = hexbin_data_object, column_start_text = column_start_text)

  }

  ## high-D model only contains the bins in 2D
  df_bin <- df_bin |>
    dplyr::filter(hb_id %in% df_bin_centroids$hexID)

  return(list(df_bin = df_bin, df_bin_centroids = df_bin_centroids))

}

find_benchmark_value <- function(.data, distance_col) {

  .data <- .data |>
    dplyr::mutate(dplyr::across({
      {
        distance_col
      }
    }, \(x) round(x, 3)))


  sorted_distance_df <- .data |>
    dplyr::arrange({
      {
        distance_col
      }
    })  ## Sort the distances

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

compute_mean_density_hex <- function(df_bin_centroids, num_bins_x) {

  # To store mean densities of hexagons
  mean_density_vec <- c()

  for (i in 1:length(df_bin_centroids$hexID)) {

    ## Identify neighbors of a specific hex bin
    neighbor_df <- df_bin_centroids |>
      dplyr::filter((hexID == (df_bin_centroids$hexID[i] + 1)) | (hexID == (df_bin_centroids$hexID[i] - 1)) |
                      (hexID == (df_bin_centroids$hexID[i] + (num_bins_x + 1))) |
                      (hexID == (df_bin_centroids$hexID[i] + num_bins_x)) |
                      (hexID == (df_bin_centroids$hexID[i] - (num_bins_x + 1))) |
                      (hexID == (df_bin_centroids$hexID[i] - num_bins_x)))

    mean_density <- neighbor_df |>
      dplyr::pull(std_counts) |>
      sum()/NROW(neighbor_df) ## The reason to take the mean is to check the density in a considerable amount

    mean_density_vec <- append(mean_density_vec, mean_density)

  }

  df_bin_centroids <- df_bin_centroids |>
    dplyr::mutate(mean_density = mean_density_vec)

  return(df_bin_centroids)

}

find_low_density_hexagons <- function(df_bin_centroids_all, num_bins_x, df_bin_centroids_low) {
  ## To compute mean density of hexagons
  df_bin_centroids <- compute_mean_density_hex(df_bin_centroids_all, num_bins_x)
  mean_density_vec <- df_bin_centroids$mean_density

  df_bin_centroids_low <- df_bin_centroids |>
    dplyr::filter(hexID %in% df_bin_centroids_low$hexID)

  ## Take first quartile
  benchmark_mean_dens_rm_hex <- stats::quantile(mean_density_vec, probs = c(0,0.25,0.5,0.75,1))[2]

  remove_bins <- c()

  ## Check only already identified low-density hexagons
  for (i in 1:length(df_bin_centroids_low$hexID)) {

    df_bin_centroids_coordinates_spec_bin <- df_bin_centroids_low |>
      dplyr::filter(hexID == df_bin_centroids_low$hexID[i])

    bin_ID <- df_bin_centroids_coordinates_spec_bin |>
      dplyr::pull(hexID)


    if(df_bin_centroids_coordinates_spec_bin$mean_density < benchmark_mean_dens_rm_hex){
      remove_bins <- append(remove_bins, bin_ID)
    }
  }

  return(remove_bins)
}

extract_coord_of_shifted_hex_grid <- function(nldr_data_with_hb_id, num_bins_x, hex_full_count_df, shift = NA) {

  if (is.na(shift)) {
    cell_diameter <- sqrt(2 * 1 / sqrt(3))
    shift <- cell_diameter/2

  }

  ## Filter centroids with their hexIDs
  hexbin_coord_all <- hex_full_count_df |>
    dplyr::select(c_x, c_y, hexID) |>
    dplyr::distinct()

  hexbin_coord_all_new <- hexbin_coord_all |>
    dplyr::mutate(c_x = c_x - shift,
                  c_y = c_y - shift) |>
    dplyr::rename(c("x" = "c_x",
                    "y" = "c_y"))

  ## Generate all coordinates of hexagons
  hex_grid_new <- full_hex_grid(hexbin_coord_all_new)

  hexbin_coord_all_new <- hexbin_coord_all_new |>
    dplyr::rename(c("c_x" = "x",
                    "c_y" = "y"))

  ## Map the polygon ID to the hexagon coordinates
  full_grid_with_polygon_id_df <- map_polygon_id(hexbin_coord_all_new, hex_grid_new)

  full_grid_with_hexbin_id_rep <- full_grid_with_polygon_id_df |>
    dplyr::slice(rep(1:dplyr::n(), each = 6)) |>
    dplyr::arrange(polygon_id)

  ## Generate the dataset with polygon, and hexagon bin centroid coordinates
  hex_full_count_df_new <- dplyr::bind_cols(hex_grid_new, full_grid_with_hexbin_id_rep)

  ## Datafarme to store new hexIDs
  nldr_df_with_new_hexID <- data.frame(matrix(ncol = 0, nrow = 0))

  for (i in 1:NROW(nldr_data_with_hb_id)) {

    ## Select the nldr point
    nldr_data_with_hb_id_spec <- nldr_data_with_hb_id |>
      dplyr::filter(dplyr::row_number() == i)

    ## Find nearest hexIDs
    df_bin_centroids_coordinates_spec_bin_near1 <- hexbin_coord_all_new |>
      dplyr::filter((hexID == nldr_data_with_hb_id_spec$hb_id[1]) |(hexID == (nldr_data_with_hb_id_spec$hb_id[1] + (num_bins_x + 1))) | (hexID == (nldr_data_with_hb_id_spec$hb_id[1] + num_bins_x)) | (hexID == (nldr_data_with_hb_id_spec$hb_id[1] - (num_bins_x + 1))) | (hexID == (nldr_data_with_hb_id_spec$hb_id[1] - num_bins_x)))

    nldr_data_with_hb_id_spec <- nldr_data_with_hb_id_spec |>
      dplyr::select(-ID) |>
      dplyr::rename("x" = names(nldr_data_with_hb_id_spec)[1],
                    "y" = names(nldr_data_with_hb_id_spec)[2])

    df_bin_centroids_coordinates_spec_bin_near1 <- df_bin_centroids_coordinates_spec_bin_near1 |>
      dplyr::rename("x" = "c_x",
                    "y" = "c_y",
                    "hb_id" = "hexID")

    near_df_1 <- dplyr::bind_rows(nldr_data_with_hb_id_spec, df_bin_centroids_coordinates_spec_bin_near1)

    ## Compute the distance from selected point to neighbouring centroids
    near_df_1$distance <- lapply(seq(nrow(near_df_1)), function(x) {
      start <- unlist(near_df_1[1, c("x","y")])
      end <- unlist(near_df_1[x, c("x","y")])
      sqrt(sum((start - end)^2))})

    near_df_1$distance <- unlist(near_df_1$distance)

    near_df_1 <- near_df_1 |>
      dplyr::filter(dplyr::row_number() != 1) |>
      dplyr::arrange(distance)

    ## Select the most nearest centroid and assign the hexID of that centroid
    nldr_data_with_hb_id_spec <- nldr_data_with_hb_id_spec |>
      dplyr::select(-hb_id) |>
      dplyr::mutate(hb_id = near_df_1$hb_id[1])

    nldr_df_with_new_hexID <- dplyr::bind_rows(nldr_df_with_new_hexID, nldr_data_with_hb_id_spec)

  }


  ## Find counts within each hexagon
  hb_id_with_counts <- nldr_df_with_new_hexID |>
    dplyr::count(hb_id) |>
    dplyr::mutate(counts = n,
                  std_counts = n/max(n)) |>
    dplyr::select(-n)

  hex_full_count_df_new <- dplyr::left_join(hex_full_count_df_new, hb_id_with_counts,
                                            by = c("hexID" = "hb_id"))

  return(hex_full_count_df_new)

}

compute_aic <- function(p, total, num_bins, num_obs) {
  mse <- mean(total) / p
  aic <- 2*num_bins*p + num_obs*p*log(mse)
  return(aic)
}

predict_hex_id <- function(df_bin_centroids, nldr_df_test, x = "UMAP1", y = "UMAP2") {

  df_bin_centroids <- df_bin_centroids |>
    dplyr::select(x, y, hexID)

  pred_hb_id <- class::knn(df_bin_centroids |> dplyr::select(-hexID),
                           nldr_df_test |> dplyr::select(!!! rlang::syms(c(x, y))),
                           cl = df_bin_centroids$hexID)

  pred_data <- nldr_df_test |>
    dplyr::mutate(pred_hb_id = as.numeric(as.character(pred_hb_id)))

  return(pred_data)

}

generate_eval_df <- function(data, prediction_df, df_bin_centroids, df_bin,
                             num_bins, col_start = "x") {

  ## Generate all possible bin centroids in the full grid
  full_centroid_df <- generate_full_grid_centroids(df_bin_centroids)

  df_bin_centroids_filtered <- df_bin_centroids |>
    dplyr::select(hexID, x, y)

  ## To map centroid coordinates to predicted hexID
  prediction_df <- dplyr::inner_join(prediction_df, df_bin_centroids_filtered,
                                     by = c("pred_hb_id" = "hexID"))

  df_bin_train <- df_bin
  names(df_bin_train)[-1] <- paste0("avg_", names(df_bin_train)[-1])

  prediction_df <- prediction_df |>
    dplyr::left_join(df_bin_train, by = c("pred_hb_id" = "hb_id")) ## Map high-D averaged/weighted mean coordinates

  prediction_df <- prediction_df |>
    dplyr::left_join(data, by = c("ID" = "ID")) ## Map high-D data

  for (i in 1:(NCOL(df_bin_train) - 1)) {

    prediction_df[ , paste0("error_square_", col_start, i)] <- (prediction_df[ , paste0(col_start, i)] - prediction_df[ , paste0("avg_", col_start, i)])^2

  }

  prediction_df <- prediction_df |>
    dplyr::mutate(total = rowSums(dplyr::pick(tidyselect::starts_with(paste0("error_square_", col_start)))))


  #number_of_bins: Total number of bins with empty bins
  eval_df <- tibble::tibble(number_of_bins = NROW(full_centroid_df),
                            number_of_observations = NROW(prediction_df),
                            total_error = compute_aic((NCOL(df_bin) - 1), prediction_df$total, NROW(df_bin_centroids), NROW(prediction_df)),
                            total_mse = mean(prediction_df$total))

  return(eval_df)

}

predict_2d_embeddings <- function(test_data, df_bin_centroids, df_bin, type_NLDR = "UMAP") {

  columns_df <- c(paste0("pred_", type_NLDR, "_", 1:2), "ID", "pred_hb_id")
  vec <- stats::setNames(rep("", length(columns_df)), columns_df)  ## Define column names

  predict_coord_test <- dplyr::bind_rows(vec)[0, ]
  predict_coord_test <- predict_coord_test |>
    dplyr::mutate_if(is.character, as.numeric)

  for (i in 1:NROW(test_data)) {

    ### Filter the new data point
    test_data_point <- test_data |>
      dplyr::filter(dplyr::row_number() == i)

    ## Obtain centroid coordinates in high-D
    centroid_coord_high_D <- df_bin |>
      dplyr::select(-hb_id)

    ## Compute the distance between test point and the centroid points in high-D
    d <- stats::dist(dplyr::bind_rows(test_data_point |> dplyr::select(-ID), centroid_coord_high_D)) |> as.matrix()

    ## Obtain the distances
    distance_vec <- d[2:dim(d)[1], 1] |> as.vector()

    ## Add the distance vec as a column in high-D centroid coordinate data set
    centroid_coord_high_D <- centroid_coord_high_D |>
      dplyr::mutate(distance = distance_vec) |>
      dplyr::mutate(hb_id = df_bin$hb_id)

    ## Sort by distance and obtain the centroid which is nearest
    predict_centroid_coord_high_D <- centroid_coord_high_D |>
      dplyr::arrange(distance) |>
      dplyr::filter(dplyr::row_number() == 1)

    ## Rename columns
    #names(predict_centroid_coord_high_D)[1:(NCOL(test_data_point) - 1)] <- paste0("C_", names(predict_centroid_coord_high_D)[1:(NCOL(test_data_point) - 1)])

    ## Obtain 2D coordinate of the nearest high-D centroid
    predict_centroid_coord_2D <- df_bin_centroids |>
      dplyr::filter(hexID %in% predict_centroid_coord_high_D$hb_id) |>
      dplyr::select(x, y) |>
      dplyr::mutate(ID = test_data_point$ID,
                    pred_hb_id = predict_centroid_coord_high_D$hb_id)

    ## Rename columns
    names(predict_centroid_coord_2D) <- c(paste0("pred_", type_NLDR, "_", 1:2), "ID", "pred_hb_id")

    ## Combine high-D and 2D coordinate
    #predict_centroid_coord_all <- dplyr::bind_cols(test_data_point, predict_centroid_coord_high_D, predict_centroid_coord_2D)

    ## Combine all
    predict_coord_test <- dplyr::bind_rows(predict_coord_test, predict_centroid_coord_2D)


  }

  return(predict_coord_test)


}

