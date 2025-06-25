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

df <- dplyr::bind_rows(true_model, training_data)

langevitour(df |> dplyr::select(-type), lineFrom = append(point_connect_df$from, point_connect_df_dir_2$from),
            lineTo = append(point_connect_df$to, point_connect_df_dir_2$to),
            group = df$type, pointSize = append(rep(2, NROW(true_model)), rep(0.8, NROW(training_data))), levelColors = c("#6a3d9a", "#969696"),
            lineColors = rep("#969696", length(append(point_connect_df$to, point_connect_df_dir_2$from))))


######################with high-D model
library(ggplot2)
library(ggbeeswarm)
source("quollr_code.R", local = TRUE)
## Import data
training_data <- read_rds("data/s_curve/s_curve_training.rds")
UMAP_s_curve <- read_rds("data/s_curve/s_curve_umap.rds")

## UMAP

num_bins_umap_s_curve <- 8
shape_val_umap_s_curve <- calculate_effective_shape_value(.data = UMAP_s_curve,
                                                          x = UMAP1, y = UMAP2) ## 1.259938
## To extract bin centroids
hexbin_data_object_umap_s_curve <- extract_hexbin_centroids(nldr_df = UMAP_s_curve, num_bins = num_bins_umap_s_curve, shape_val = shape_val_umap_s_curve, x = UMAP1, y = UMAP2)

df_bin_centroids_umap_s_curve <- hexbin_data_object_umap_s_curve$hexdf_data

UMAP_data_with_hb_id_s_curve <- UMAP_s_curve |>
  dplyr::mutate(hb_id = hexbin_data_object_umap_s_curve$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_umap_s_curve <- dplyr::bind_cols(training_data |> dplyr::select(-ID), UMAP_data_with_hb_id_s_curve)

## Averaged on high-D
df_bin_umap_s_curve <- avg_highD_data(.data = df_all_umap_s_curve)

## Triangulate bin centroids
tr1_object_umap_s_curve <- triangulate_bin_centroids(df_bin_centroids_umap_s_curve, x, y)
tr_from_to_df_umap_s_curve <- generate_edge_info(triangular_object = tr1_object_umap_s_curve)

# ggplot(df_bin_centroids_umap_s_curve, aes(x = x, y = y)) +
#   geom_segment(data = tr_from_to_df_umap_s_curve, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
#   geom_point(size = 2, colour = "#33a02c") +
#   coord_equal()


## Compute 2D distances
distance_umap_s_curve <- cal_2d_dist(.data = tr_from_to_df_umap_s_curve)

plot_dist <- function(distance_df){
  distance_df$group <- "1"
  dist_plot <- ggplot(distance_df, aes(x = group, y = distance)) +
    geom_quasirandom()+
    ylim(0, max(unlist(distance_df$distance))+ 0.5) + coord_flip()
  return(dist_plot)
}

plot_dist(distance_umap_s_curve)

## To find the benchmark value
#benchmark_umap_s_curve <- find_benchmark_value(.data = distance_umap_s_curve, distance_col = distance)
benchmark_umap_s_curve <- 1.346694

df_b <- df_bin_umap_s_curve |>
  dplyr::filter(hb_id %in% df_bin_centroids_umap_s_curve$hexID) |>
  dplyr::select(-hb_id) |>
  dplyr::mutate(type = "model")

distance_df_small_edges <- distance_umap_s_curve %>%
  dplyr::filter(distance < benchmark_umap_s_curve)

distance_df_small_edges$from <- distance_df_small_edges$from + max(point_connect_df_dir_2$from)
distance_df_small_edges$to <- distance_df_small_edges$to + max(point_connect_df_dir_2$to)

point_connect_all <- dplyr::bind_rows(point_connect_df, point_connect_df_dir_2, distance_df_small_edges |> dplyr::select(from, to))

model_df <- dplyr::bind_rows(true_model, df_b)

df <- dplyr::bind_rows(model_df, training_data |> select(-ID) |> mutate(type = "data"))

langevitour(df |> dplyr::select(-type), lineFrom = point_connect_all$from,
            lineTo = point_connect_all$to,
            group = df$type, pointSize = append(rep(1, NROW(true_model)), rep(3, NROW(df_b) + NROW(training_data))), levelColors = c("#6a3d9a", "#33a02c", "#d95f02"),
            lineColors = append(rep("#e41a1c", length(append(point_connect_df$from, point_connect_df_dir_2$from))), rep("black", length(distance_df_small_edges$from))))



##########################

# langevitour::langevitour(training_data |> dplyr::select(-ID) |> dplyr::select(x1,x3))
#
# langevitour::langevitour(training_data |> dplyr::select(-ID) |> dplyr::select(x1,x2))
#
# langevitour::langevitour(training_data |> dplyr::select(-ID) |> dplyr::select(x2,x3))
#
#
# min(training_data$x1)
# max(training_data$x1)
#
# seq_points_x1 <- seq(min(training_data$x1), max(training_data$x1), by = 0.3)
#
# min(training_data$x2)
# max(training_data$x2)
#
# seq_points_x2 <- seq(min(training_data$x2), max(training_data$x2), by = 0.3)
#
# x1_x2_wire_points <- expand.grid(seq_points_x1, seq_points_x2)
# names(x1_x2_wire_points) <- c("x1", "x2")
#
#
# corner_points <- tibble::tibble(x1 = c(min(training_data$x1), min(training_data$x1), max(training_data$x1), max(training_data$x1)),
#                                 x2 = c(min(training_data$x2), max(training_data$x2),
#                                       min(training_data$x2), max(training_data$x2)))
#
# langevitour::langevitour(x1_x2_wire_points)
#
# seq_points_x3 <- seq(min(training_data$x3), max(training_data$x3), by = 0.3)
#
# x1_x3_wire_points <- expand.grid(seq_points_x1, seq_points_x3)
# names(x1_x3_wire_points) <- c("x1", "x3")
# langevitour::langevitour(x1_x3_wire_points)
#
# a <- training_data |>
#   dplyr::filter(x1 %in% quantile(training_data$x1)) |>
#   dplyr::pull(x1)
#
# b <- training_data |>
#   dplyr::filter(x3 %in% quantile(training_data$x3)) |>
#   dplyr::pull(x3)
#
# x1_x3_wire_points <- expand.grid(a, b)
# names(x1_x3_wire_points) <- c("x1", "x3")
#
# langevitour::langevitour(x1_x3_wire_points)
#
# bid <- expand.grid(x1_x2_wire_points$x1, x1_x2_wire_points$x2, training_data$x3 |> quantile())
# langevitour::langevitour(bid)
#
# #####################
# training_data_n <- training_data |> dplyr::mutate(asin_x1 = asin(x1))
#
# min_max_seq_points_asin_x1 <- append(min(training_data_n$asin_x1), max(training_data_n$asin_x1))
#
# seq_points_asin_x1 <- sample(training_data_n$asin_x1[!(training_data_n$asin_x1 %in% min_max_seq_points_asin_x1)], 20)
#
# filtered_training_data_n <- training_data_n |>
#   dplyr::filter(asin_x1 %in% seq_points_asin_x1)
#
# tibble::tibble(x1 = filtered_training_data_n$x1, )
#
# min_max_seq_points_asin_x3 <- append(min(training_data_n$asin_x1), max(training_data_n$asin_x1))
# seq_points_asin_x3 <- sample(training_data_n$asin_x1[!(training_data_n$asin_x1 %in% min_max_seq_points_asin_x3)], 20)
#
# filtered_training_data_n <- training_data_n |>
#   dplyr::filter(asin_x1 %in% seq_points_asin_x3)
#
# ########
# tt <- 3 * pi * stats::runif(n = 100, min = -0.5, max = 0.5)
#
# tt <- training_data_n$asin_x1 |> unique() |> sort() |> sample(20)
# x <- sin(tt)
# y <- 2.0 * stats::runif(n = 20)
# z <- sign(tt) * (cos(tt) - 1)
# X <- cbind(x, y, z)
# langevitour(X)
#
# ##### Tried with interpolation
# training_data_sample <- read_rds("data/s_curve/s_curve_training.rds") |>
#   sample_n(size = 100)
#
# langevitour::langevitour(training_data_sample |> dplyr::select(-ID))
#
# # Apply approx function
# data_approx1 <- approx(training_data_sample$x1, training_data_sample$x3)
# data_approx1
#
# # Draw output of approx function
# plot(data_approx1$x,
#      data_approx1$y)
# points(training_data_sample$x1, training_data_sample$x2,
#        col = "red",
#        pch = 16)
#
# training_data <- read_rds("data/s_curve/s_curve_training.rds")
#
# seq_points_x1 <- seq(min(training_data$x1), max(training_data$x1), by = 0.1)
#
# # Perform spline interpolation
# interp_y <- spline(training_data$x1, training_data$x3, xout = seq_points_x1)$y
#
# ggplot(data = tibble::tibble(seq_points_x1, interp_y), aes(x = seq_points_x1, y = interp_y)) + geom_point()
# langevitour::langevitour(tibble::tibble(seq_points_x1, interp_y))

# # Install and load the deldir package
# #install.packages("deldir")
# library(deldir)
#
# # Generate some sample points
# set.seed(123)
# points <- df_bin_centroids_umap_s_curve |> dplyr::select(x, y)
#
# # Compute the Delaunay triangulation and Voronoi diagram
# vd <- deldir(points[,1], points[,2])
#
# # Plot the Voronoi diagram
# plot(vd)
