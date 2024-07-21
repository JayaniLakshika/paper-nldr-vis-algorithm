## Import necessary libraries
library(quollr)
library(dplyr)
library(reader)
library(langevitour)

training_data_scurve <- read_rds("data/s_curve/scurve_500.rds")

tsne_scurve <- read_rds(file = "data/s_curve/new/s_curve_tsne_22.rds")

scurve_scaled_obj <- gen_scaled_data(
  data = tsne_scurve)

tsne_scurve_scaled <- scurve_scaled_obj$scaled_nldr
lim1 <- scurve_scaled_obj$lim1
lim2 <- scurve_scaled_obj$lim2
r2 <- diff(lim2)/diff(lim1)

## Compute hexbin parameters
num_bins_x_scurve <- 12

scurve_model <- fit_highd_model(
  training_data = training_data_scurve,
  emb_df = tsne_scurve_scaled,
  bin1 = num_bins_x_scurve,
  r2 = r2,
  q = 0.12,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = TRUE,
  col_start_highd = "x"
)

df_bin_centroids_scurve <- scurve_model$df_bin_centroids
df_bin_scurve <- scurve_model$df_bin

## Triangulate bin centroids
tr1_object_scurve <- tri_bin_centroids(
  df_bin_centroids_scurve, x = "c_x", y = "c_y")
tr_from_to_df_scurve <- gen_edges(
  tri_object = tr1_object_scurve)

## Compute 2D distances
distance_scurve <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_scurve,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_scurve <- find_lg_benchmark(
  distance_edges = distance_scurve,
  distance_col = "distance")

trimesh_removed_scurve <- vis_rmlg_mesh(
  distance_edges = distance_scurve,
  benchmark_value = benchmark_scurve,
  tr_coord_df = tr_from_to_df_scurve,
  distance_col = "distance")

## Hexagonal binning to have regular hexagons
hb_obj_scurve <- hex_binning(
  data = tsne_scurve_scaled,
  bin1 = num_bins_x_scurve,
  r2 = r2,
  q = 0.12)

# ## Data set with all possible centroids in the hexagonal grid
# all_centroids_df <- hb_obj_scurve$centroids
# glimpse(all_centroids_df)
#
# ## Generate all coordinates of hexagons
# hex_grid <- hb_obj_scurve$hex_poly
# glimpse(hex_grid)
#
# ## To obtain the standardise counts within hexbins
# counts_df <- hb_obj_scurve$std_cts
# df_bin_centroids <- extract_hexbin_centroids(centroids_df = all_centroids_df,
#                                              counts_df = counts_df) |>
#   filter(drop_empty == FALSE)
#
# hex_grid_with_counts <- left_join(hex_grid, counts_df, by = c("hex_poly_id" = "hb_id"))
#
# ggplot(data = hex_grid_with_counts, aes(x = x, y = y)) +
#   geom_polygon(color = "black", aes(group = hex_poly_id, fill = std_counts)) +
#   geom_text(data = all_centroids_df, aes(x = c_x, y = c_y, label = hexID)) +
#   scale_fill_viridis_c(direction = -1, na.value = "#ffffff") +
#   coord_fixed()

tsne_data_with_hb_id <- hb_obj_scurve$data_hb_id

df_all_scurve <- dplyr::bind_cols(training_data_scurve |> dplyr::select(-ID),
                                  tsne_data_with_hb_id)

### Define type column
df <- df_all_scurve |>
  dplyr::select(tidyselect::starts_with("x")) |>
  dplyr::mutate(type = "data") ## original dataset

df_b <- df_bin_scurve |>
  dplyr::filter(hb_id %in% df_bin_centroids_scurve$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_scurve$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

df_exe <- dplyr::bind_rows(df_b, df)

## Set the maximum difference as the criteria
distance_df_small_edges <- distance_scurve |>
  dplyr::filter(distance < benchmark_scurve)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to,
                         group = df_exe$type, pointSize = append(rep(0, NROW(df_b)), rep(0.5, NROW(df))),
                         levelColors = c("#6a3d9a", "#33a02c"))

##############################

### With true model

training_data <- read_rds("data/s_curve/s_curve_training.rds") |>
  dplyr::select(-ID) |>
  dplyr::mutate(type = "data")

corner_points <- tibble::tibble(x1 = c(min(training_data$x1), min(training_data$x1), ## added
                                       max(training_data$x1), max(training_data$x1), ## added
                                       mean(training_data$x1), mean(training_data$x1),
                                       mean(training_data$x1), mean(training_data$x1),
                                       mean(training_data$x1), mean(training_data$x1),
                                       sort(training_data$x1)[210], sort(training_data$x1)[210],
                                       sort(training_data$x1)[537], sort(training_data$x1)[537],
                                       sort(training_data$x1)[665], sort(training_data$x1)[665],
                                       sort(training_data$x1)[98], sort(training_data$x1)[98],
                                       sort(training_data$x1)[280], sort(training_data$x1)[280],
                                       sort(training_data$x1)[665], sort(training_data$x1)[537],
                                       sort(training_data$x1)[665], sort(training_data$x1)[537],
                                       sort(training_data$x1)[665], sort(training_data$x1)[537]),
                                x2 = c(min(training_data$x2), max(training_data$x2),
                                       min(training_data$x2), max(training_data$x2),
                                       min(training_data$x2), max(training_data$x2),
                                       min(training_data$x2), max(training_data$x2),
                                       min(training_data$x2), max(training_data$x2),
                                       min(training_data$x2), max(training_data$x2),
                                       min(training_data$x2), max(training_data$x2),
                                       min(training_data$x2), max(training_data$x2),
                                       min(training_data$x2), max(training_data$x2),
                                       min(training_data$x2), max(training_data$x2),
                                       mean(training_data$x2), mean(training_data$x2),
                                       min(training_data$x2), min(training_data$x2),
                                       max(training_data$x2), max(training_data$x2)),
                                x3 = c(training_data |> dplyr::filter(x1 == min(x1)) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == min(x1)) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == max(x1)) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == max(x1)) |> dplyr::pull(x3),
                                       min(training_data$x3),
                                       min(training_data$x3),
                                       max(training_data$x3),
                                       max(training_data$x3),
                                       0, 0,
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[210]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[210]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[537]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[537]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[665]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[665]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[98]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[98]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[280]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[280]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[665]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[537]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[665]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[537]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[665]) |> dplyr::pull(x3),
                                       training_data |> dplyr::filter(x1 == sort(training_data$x1)[537]) |> dplyr::pull(x3)
                                ))


point_connect_df <- tibble::tibble(from = c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25),
                                   to = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26))

langevitour(corner_points, lineFrom = point_connect_df$from, lineTo = point_connect_df$to)

true_model <- corner_points |>
  dplyr::mutate(type = "true model")

# true_model <- corner_points |>
#   dplyr::mutate(type = as.character(c(1:NROW(corner_points))))

df <- dplyr::bind_rows(true_model, training_data |> dplyr::select(x1:x3, type))

langevitour(df |> dplyr::select(-type), lineFrom = point_connect_df$from,
            lineTo = point_connect_df$to,
            group = df$type, pointSize = 3, levelColors = c("#6a3d9a", "#d95f02"),
            lineColors = rep("#e41a1c", NROW(point_connect_df)))

######### Develop the automate function (CORRECT)

## To create the true model for S-curve containing three parts
library(readr)
library(dplyr)
library(langevitour)

training_data <- read_rds("data/s_curve/s_curve_training.rds") |>
  dplyr::select(-ID) |>
  dplyr::mutate(type = "data")

## 1) For min and max values of x2 we need to find several options with x1 and x3
### (i) with min and max values of x1 and their corresponding x3 values
corner_points_min_max <- tibble::tibble(x1 = c(min(training_data$x1), min(training_data$x1),
                                               max(training_data$x1), max(training_data$x1)),
                                        x2 = c(min(training_data$x2), max(training_data$x2),
                                               min(training_data$x2), max(training_data$x2)),
                                        x3 = c(training_data |> dplyr::filter(x1 == min(x1)) |> dplyr::pull(x3),
                                               training_data |> dplyr::filter(x1 == min(x1)) |> dplyr::pull(x3),
                                               training_data |> dplyr::filter(x1 == max(x1)) |> dplyr::pull(x3),
                                               training_data |> dplyr::filter(x1 == max(x1)) |> dplyr::pull(x3)))

point_connect_df <- tibble::tibble(from = seq(1, NROW(corner_points_min_max), by = 2),
                                   to = seq(2, NROW(corner_points_min_max), by = 2))

langevitour(corner_points_min_max, lineFrom = point_connect_df$from, lineTo = point_connect_df$to)

true_model <- corner_points_min_max |>
  dplyr::mutate(type = "true model")

df <- dplyr::bind_rows(true_model, training_data |> dplyr::select(x1:x3, type))

langevitour(df |> dplyr::select(-type), lineFrom = point_connect_df$from,
            lineTo = point_connect_df$to,
            group = df$type, pointSize = 3, levelColors = c("#6a3d9a", "#d95f02"),
            lineColors = rep("#e41a1c", NROW(point_connect_df)))

### (ii) with mean values of x1 and their corresponding x3 values

corner_points_mean <- tibble::tibble(x1 = c(mean(training_data$x1), mean(training_data$x1),
                                            mean(training_data$x1), mean(training_data$x1),
                                            mean(training_data$x1), mean(training_data$x1)),
                                     x2 = c(min(training_data$x2), max(training_data$x2),
                                            min(training_data$x2), max(training_data$x2),
                                            min(training_data$x2), max(training_data$x2)),
                                     x3 = c(min(training_data$x3),
                                            min(training_data$x3),
                                            max(training_data$x3),
                                            max(training_data$x3),
                                            0, 0))

point_connect_df <- tibble::tibble(from = seq(1, NROW(corner_points_mean), by = 2),
                                   to = seq(2, NROW(corner_points_mean), by = 2))

langevitour(corner_points_mean, lineFrom = point_connect_df$from, lineTo = point_connect_df$to)

true_model <- corner_points_mean |>
  dplyr::mutate(type = "true model")

df <- dplyr::bind_rows(true_model, training_data |> dplyr::select(x1:x3, type))

langevitour(df |> dplyr::select(-type), lineFrom = point_connect_df$from,
            lineTo = point_connect_df$to,
            group = df$type, pointSize = 3, levelColors = c("#6a3d9a", "#d95f02"),
            lineColors = rep("#e41a1c", NROW(point_connect_df)))



## 2) For quantile values of x1 we need to find several options with x1 and x3

tolerance <- 0.001 # Tolerance level
values <- training_data$x1

## To store approx quantile values
quantile_df_all <- data.frame(matrix(nrow = 0, ncol = 0))

for (i in 1:length(quantile(values, probs = seq(.1, .9, by = .1)))) {

  specific_value <- quantile(values, probs = seq(.1, .9, by = .1))[i] # Specific value to compare against

  # Find approximately equal values
  approximately_equal <- values >= (specific_value - tolerance) &
    values <= (specific_value + tolerance)

  quantile_df_spec <- training_data |> dplyr::filter(x1 == values[approximately_equal][1])

  quantile_df_all <- dplyr::bind_rows(quantile_df_all, quantile_df_spec)

}

tolerance <- 0.003 # Tolerance level
values <- training_data$x1[training_data$x1 > 0]

## To store approx quantile values
quantile_df_greater_than_0 <- data.frame(matrix(nrow = 0, ncol = 0))

for (i in 1:length(quantile(values, probs = seq(.1, .9, by = .1)))) {

  specific_value <- quantile(values, probs = seq(.1, .9, by = .1))[i] # Specific value to compare against

  # Find approximately equal values
  approximately_equal <- values >= (specific_value - tolerance) &
    values <= (specific_value + tolerance)

  quantile_df_spec <- training_data |> dplyr::filter(x1 == values[approximately_equal][1])

  quantile_df_greater_than_0 <- dplyr::bind_rows(quantile_df_greater_than_0, quantile_df_spec)

}

tolerance <- 0.008
values <- training_data$x1[training_data$x1 < 0]

## To store approx quantile values
quantile_df_less_than_0 <- data.frame(matrix(nrow = 0, ncol = 0))

for (i in 1:length(quantile(values, probs = seq(.1, .9, by = .1)))) {

  specific_value <- quantile(values, probs = seq(.1, .9, by = .1))[i] # Specific value to compare against

  # Find approximately equal values
  approximately_equal <- values >= (specific_value - tolerance) &
    values <= (specific_value + tolerance)

  quantile_df_spec <- training_data |> dplyr::filter(x1 == values[approximately_equal][1])

  quantile_df_less_than_0 <- dplyr::bind_rows(quantile_df_less_than_0, quantile_df_spec)

}

quantile_df <- dplyr::bind_rows(quantile_df_less_than_0, quantile_df_greater_than_0)

corner_points_quantiles <- tibble::tibble(x1 = rep(quantile_df$x1, each = 2),
                                          x2 = rep(c(min(training_data$x2), max(training_data$x2)), NROW(quantile_df)),
                                          x3 = rep(quantile_df$x3, each = 2))

point_connect_df <- tibble::tibble(from = seq(1, NROW(corner_points_quantiles), by = 2),
                                   to = seq(2, NROW(corner_points_quantiles), by = 2))

langevitour(corner_points_quantiles, lineFrom = point_connect_df$from, lineTo = point_connect_df$to)

true_model <- corner_points_quantiles |>
  dplyr::mutate(type = "true model")

df <- dplyr::bind_rows(true_model, training_data |> dplyr::select(x1:x3, type))

langevitour(df |> dplyr::select(-type), lineFrom = point_connect_df$from,
            lineTo = point_connect_df$to,
            group = df$type, pointSize = 3, levelColors = c("#6a3d9a", "#d95f02"),
            lineColors = rep("#e41a1c", NROW(point_connect_df)))

# ## 3) For quantile values of x2 we need to find several options with x1 and x3
#
# tolerance <- 0.003 # Tolerance level
# values <- training_data$x2
#
# ## To store approx quantile values
# quantile_df_all_x2 <- data.frame(matrix(nrow = 0, ncol = 0))
#
# for (i in 1:length(quantile(values, probs = c(0,0.25,0.5,0.75,1)))) {
#
#   specific_value <- quantile(values, probs = c(0,0.25,0.5,0.75,1))[i] # Specific value to compare against
#
#   # Find approximately equal values
#   approximately_equal <- values >= (specific_value - tolerance) &
#     values <= (specific_value + tolerance)
#
#   quantile_df_spec <- training_data |> dplyr::filter(x2 == values[approximately_equal][1])
#
#   quantile_df_all_x2 <- dplyr::bind_rows(quantile_df_all_x2, quantile_df_spec)
#
# }
#
# corner_points_quantiles_x2 <- tibble::tibble(x1 = rep(quantile_df_all_x2$x1, each = 5),
#                                           x2 = rep(quantile_df_all_x2$x2, length(quantile_df_all_x2$x1)),
#                                           x3 = rep(quantile_df_all_x2$x3, each = 5)) |>
#   dplyr::arrange(x2)
#
# point_connect_df <- tibble::tibble(from = seq(2, (NROW(corner_points_quantiles_x2) - 2), by = 1)[-NROW(corner_points_quantiles_x2)],
#                                    to = seq(3, (NROW(corner_points_quantiles_x2) - 1), by = 1)[-NROW(corner_points_quantiles_x2)])
#
# langevitour(corner_points_quantiles_x2, lineFrom = c(1), lineTo = c(3))
#
# true_model <- corner_points_quantiles_x2 |>
#   dplyr::mutate(type = "true model")
#
# df <- dplyr::bind_rows(true_model, training_data |> dplyr::select(x1:x3, type))
#
# langevitour(df |> dplyr::select(-type), lineFrom = point_connect_df$from,
#             lineTo = point_connect_df$to,
#             group = df$type, pointSize = 3, levelColors = c("#6a3d9a", "#d95f02"),
#             lineColors = rep("#e41a1c", NROW(point_connect_df)))
#


## Final (TRUE MODEL) - one direction

corner_points_df <- dplyr::bind_rows(corner_points_min_max, corner_points_quantiles, corner_points_mean)

point_connect_df <- tibble::tibble(from = seq(1, NROW(corner_points_df), by = 2),
                                   to = seq(2, NROW(corner_points_df), by = 2))

langevitour(corner_points_df, lineFrom = point_connect_df$from, lineTo = point_connect_df$to)

true_model <- corner_points_df |>
  dplyr::mutate(type = "true model")

df <- dplyr::bind_rows(true_model, training_data |> dplyr::select(x1:x3, type))

langevitour(df |> dplyr::select(-type), lineFrom = point_connect_df$from,
            lineTo = point_connect_df$to,
            group = df$type, pointSize = 3, levelColors = c("#6a3d9a", "#d95f02"),
            lineColors = rep("#e41a1c", NROW(point_connect_df)))

######

distinct_corner_points_df <- corner_points_df |>
  dplyr::select(x1, x3) |>
  dplyr::distinct()

corner_points_quantiles_x2_q1 <- tibble::tibble(x1 = distinct_corner_points_df$x1,
                                                x2 = rep(quantile(training_data$x2,
                                                                  probs = c(0,0.25,0.5,0.75,1))[1],
                                                         NROW(distinct_corner_points_df)),
                                                x3 = distinct_corner_points_df$x3) |>
  dplyr::arrange(x3)


new_point_connect_df_q1 <- tibble::tibble(from = c(47, 48, 49, 50, 47, 51, 54, 55, 56, 57, 58, 59, 60, 61, 64, 65, 66, 69, 68, 62, 63),
                                          to = c(48, 49, 50, 52, 51, 54, 55, 56, 57, 58, 59, 60, 61, 64, 65, 69, 67, 68, 67, 63, 66))


corner_points_quantiles_x2_q2 <- tibble::tibble(x1 = distinct_corner_points_df$x1,
                                                x2 = rep(quantile(training_data$x2,
                                                                  probs = c(0,0.25,0.5,0.75,1))[2], NROW(distinct_corner_points_df)),
                                                x3 = distinct_corner_points_df$x3) |>
  dplyr::arrange(x3)

new_point_connect_df_q2 <- tibble::tibble(from = NROW(new_point_connect_df_q1) + 2 + c(47, 48, 49, 50, 47, 51, 54, 55, 56, 57, 58, 59, 60, 61, 64, 65, 66, 69, 68, 62, 63),
                                          to = NROW(new_point_connect_df_q1) + 2 + c(48, 49, 50, 52, 51, 54, 55, 56, 57, 58, 59, 60, 61, 64, 65, 69, 67, 68, 67, 63, 66))


corner_points_quantiles_x2_q3 <- tibble::tibble(x1 = distinct_corner_points_df$x1,
                                                x2 = rep(quantile(training_data$x2,
                                                                  probs = c(0,0.25,0.5,0.75,1))[3], NROW(distinct_corner_points_df)),
                                                x3 = distinct_corner_points_df$x3) |>
  dplyr::arrange(x3)

new_point_connect_df_q3 <- tibble::tibble(from = NROW(new_point_connect_df_q1) + NROW(new_point_connect_df_q2) + 4 + c(47, 48, 49, 50, 47, 51, 54, 55, 56, 57, 58, 59, 60, 61, 64, 65, 66, 69, 68, 62, 63),
                                          to = NROW(new_point_connect_df_q1) + NROW(new_point_connect_df_q2) + 4 + c(48, 49, 50, 52, 51, 54, 55, 56, 57, 58, 59, 60, 61, 64, 65, 69, 67, 68, 67, 63, 66))



corner_points_quantiles_x2_q4 <- tibble::tibble(x1 = distinct_corner_points_df$x1,
                                                x2 = rep(quantile(training_data$x2,
                                                                  probs = c(0,0.25,0.5,0.75,1))[4], NROW(distinct_corner_points_df)),
                                                x3 = distinct_corner_points_df$x3) |>
  dplyr::arrange(x3)

new_point_connect_df_q4 <- tibble::tibble(from = NROW(new_point_connect_df_q1) + NROW(new_point_connect_df_q2) + NROW(new_point_connect_df_q3) + 6  + c(47, 48, 49, 50, 47, 51, 54, 55, 56, 57, 58, 59, 60, 61, 64, 65, 66, 69, 68, 62, 63),
                                          to = NROW(new_point_connect_df_q1) + NROW(new_point_connect_df_q2) + NROW(new_point_connect_df_q3) + 6 + c(48, 49, 50, 52, 51, 54, 55, 56, 57, 58, 59, 60, 61, 64, 65, 69, 67, 68, 67, 63, 66))



corner_points_quantiles_x2_q5 <- tibble::tibble(x1 = distinct_corner_points_df$x1,
                                                x2 = rep(quantile(training_data$x2,
                                                                  probs = c(0,0.25,0.5,0.75,1))[5], NROW(distinct_corner_points_df)),
                                                x3 = distinct_corner_points_df$x3) |>
  dplyr::arrange(x3)

new_point_connect_df_q5 <- tibble::tibble(from = NROW(new_point_connect_df_q1) + NROW(new_point_connect_df_q2) + NROW(new_point_connect_df_q3) + NROW(new_point_connect_df_q4) + 8 + c(47, 48, 49, 50, 47, 51, 54, 55, 56, 57, 58, 59, 60, 61, 64, 65, 66, 69, 68, 62, 63),
                                          to = NROW(new_point_connect_df_q1) + NROW(new_point_connect_df_q2) + NROW(new_point_connect_df_q3) + NROW(new_point_connect_df_q4) + 8 + c(48, 49, 50, 52, 51, 54, 55, 56, 57, 58, 59, 60, 61, 64, 65, 69, 67, 68, 67, 63, 66))



corner_points_quantiles_x2 <- dplyr::bind_rows(corner_points_quantiles_x2_q1, corner_points_quantiles_x2_q2,
                                               corner_points_quantiles_x2_q3, corner_points_quantiles_x2_q4,
                                               corner_points_quantiles_x2_q5)


point_connect_df_dir_2 <- bind_rows(new_point_connect_df_q1, new_point_connect_df_q2, new_point_connect_df_q3, new_point_connect_df_q4, new_point_connect_df_q5)
# #langevitour(corner_points_quantiles_x2, lineFrom = point_connect_df$from, lineTo = point_connect_df$to)
#
# true_model <- corner_points_quantiles_x2 |>
#   dplyr::mutate(type = "true model")
#
# df <- dplyr::bind_rows(true_model, training_data |> dplyr::select(x1:x3, type))
#
# langevitour(df |> dplyr::select(-type), lineFrom = point_connect_df$from,
#             lineTo = point_connect_df$to,
#             group = df$type, pointSize = 3, levelColors = c("#6a3d9a", "#d95f02"),
#             lineColors = rep("#e41a1c", NROW(point_connect_df)))

#######TRUE MODEL (final)##########

true_model <- bind_rows(corner_points_df, corner_points_quantiles_x2_q1, corner_points_quantiles_x2_q2,
                        corner_points_quantiles_x2_q3, corner_points_quantiles_x2_q4, corner_points_quantiles_x2_q5) |>
  dplyr::mutate(type = "true model")

# new_from_point_1 <- append(new_point_connect_df_q1$from, new_point_connect_df_q2$from)
# new_to_point_1 <- append(new_point_connect_df_q1$to, new_point_connect_df_q2$to)

df <- dplyr::bind_rows(true_model, training_data |> dplyr::select(x1:x3, type))

langevitour(df |> dplyr::select(-type), lineFrom = append(point_connect_df$from, point_connect_df_dir_2$from),
            lineTo = append(point_connect_df$to, point_connect_df_dir_2$to),
            group = df$type, pointSize = 3, levelColors = c("#6a3d9a", "#d95f02"),
            lineColors = rep("#e41a1c", length(append(point_connect_df$to, point_connect_df_dir_2$from))))


true_model <- true_model |>
  dplyr::mutate(x4 = mean(runif(NROW(true_model), -0.02, 0.02)),
                x5 = mean(runif(NROW(true_model), -0.02, 0.02)),
                x6 = mean(runif(NROW(true_model), -0.1, 0.1)),
                x7 = mean(runif(NROW(true_model), -0.01, 0.01)))

# df <- dplyr::bind_rows(true_model, training_data)
#
# langevitour(df |> dplyr::select(-type), lineFrom = append(point_connect_df$from, point_connect_df_dir_2$from),
#             lineTo = append(point_connect_df$to, point_connect_df_dir_2$to),
#             group = df$type, pointSize = 3, levelColors = c("#6a3d9a", "#d95f02"),
#             lineColors = rep("#e41a1c", length(append(point_connect_df$to, point_connect_df_dir_2$from))))

distance_df_small_edges$from <- distance_df_small_edges$from + max(point_connect_df_dir_2$from)
distance_df_small_edges$to <- distance_df_small_edges$to + max(point_connect_df_dir_2$to)

point_connect_all <- dplyr::bind_rows(point_connect_df, point_connect_df_dir_2, distance_df_small_edges |> dplyr::select(from, to))


model_df <- dplyr::bind_rows(true_model, df_b)

training_data_scurve <- read_rds("data/s_curve/scurve_500.rds") |>
  select(-ID)

df <- dplyr::bind_rows(model_df, training_data_scurve |> mutate(type = "data"))

langevitour(df |> dplyr::select(-type), lineFrom = point_connect_all$from,
            lineTo = point_connect_all$to,
            group = df$type, pointSize = append(rep(2, NROW(true_model) +  NROW(df_b)), rep(1, NROW(training_data_scurve))), levelColors = c("#6a3d9a", "#33a02c", "#969696"),
            lineColors = append(rep("#969696", length(append(point_connect_df$from, point_connect_df_dir_2$from))), rep("black", length(distance_df_small_edges$from))))


