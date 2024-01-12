library(readr)
library(dplyr)

set.seed(20240110)

source("quollr_code.R", local = TRUE)

data <- read_rds("data/s_curve/s_curve.rds")
training_data <- read_rds("data/s_curve/s_curve_training.rds")
test_data <- read_rds("data/s_curve/s_curve_test.rds")
UMAP_data <- read_rds(file = "data/s_curve/s_curve_umap.rds")


shape_value_curve <- calculate_effective_shape_value(.data = UMAP_data,
                                                     x = UMAP1, y = UMAP2)

## From AIC plot
num_bins_curve <- 8

## To extract bin centroids
hexbin_data_object_s_curve <- extract_hexbin_centroids(nldr_df = UMAP_data,
                                                      num_bins = num_bins_curve,
                                                      shape_val = shape_value_curve, x = UMAP1, y = UMAP2)

df_bin_centroids_s_curve <- hexbin_data_object_s_curve$hexdf_data

UMAP_s_curve_with_hb_id <- UMAP_data |>
  dplyr::mutate(hb_id = hexbin_data_object_s_curve$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_s_curve <- dplyr::bind_cols(training_data |> dplyr::select(-ID), UMAP_s_curve_with_hb_id)

## Averaged on high-D
df_bin_s_curve <- avg_highD_data(.data = df_all_s_curve, column_start_text = "x")

predict_2d_embeddings(test_data = test_data, df_bin_centroids = df_bin_centroids_s_curve, df_bin = df_bin_s_curve, type_NLDR = "UMAP")


# pred_df_training_object <- predict_hex_id(training_data = training_data, nldr_df = UMAP_data, nldr_df_test = UMAP_data, num_bins = num_bins_curve, shape_val = shape_value_curve, x = "UMAP1", y = "UMAP2", col_start = "x")
# pred_df_training <- pred_df_training_object$pred_data
# centroid_df_training <- pred_df_training_object$df_bin_centroids
# avg_df_training <- pred_df_training_object$df_bin


# eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, num_bins = num_bins_curve, col_start = "x")

### Predict embeddings for test data (Method 1)

# ### Filter the new data point
# test_data <- test_data |>
#   slice_sample(n = 1)
#
# ## Find the nearest centroid in high-D
#
# centroid_coord_high_D <- avg_df_training |>
#   select(-hb_id)
#
# d <- dist(bind_rows(test_data, centroid_coord_high_D)) |> as.matrix()
# distance_vec <- d[2:dim(d)[1], 1] |> as.vector()
#
# centroid_coord_high_D <- centroid_coord_high_D |>
#   dplyr::mutate(distance = distance_vec) |>
#   dplyr::mutate(hb_id = avg_df_training$hb_id)
#
# predict_centroid_coord_high_D <- centroid_coord_high_D |>
#   dplyr::arrange(distance) |>
#   dplyr::filter(dplyr::row_number() == 1)
#
# names(predict_centroid_coord_high_D)[1:(NCOL(test_data) - 1)] <- paste0("C_", names(predict_centroid_coord_high_D)[1:(NCOL(test_data) - 1)])
#
# predict_centroid_coord_2D <- centroid_df_training |>
#   dplyr::filter(hexID %in% predict_centroid_coord_high_D$hb_id) |>
#   dplyr::select(x, y)
#
# names(predict_centroid_coord_2D) <- paste0("C_", names(predict_centroid_coord_2D))
#
# predict_centroid_coord_all <- dplyr::bind_cols(predict_centroid_coord_high_D, predict_centroid_coord_2D)
#
# predict_coord_test <- dplyr::bind_rows(predict_coord_test, predict_centroid_coord_all)

#########
columns_df <- append(paste0("C_", names(test_data)[-NCOL(test_data)]), c("distance", "hb_id", "C_x", "C_y"))
vec <- stats::setNames(rep("", length(columns_df)), columns_df)  ## Define column names

predict_coord_test <- dplyr::bind_rows(vec)[0, ]
predict_coord_test <- predict_coord_test |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:NROW(test_data)) {

  ### Filter the new data point
  test_data_point <- test_data |>
    dplyr::filter(dplyr::row_number() == i)

  ## Obtain centroid coordinates in high-D
  centroid_coord_high_D <- avg_df_training |>
    select(-hb_id)

  ## Compute the distance between test point and the centroid points in high-D
  d <- dist(bind_rows(test_data_point, centroid_coord_high_D)) |> as.matrix()

  ## Obtain the distances
  distance_vec <- d[2:dim(d)[1], 1] |> as.vector()

  ## Add the distance vec as a column in high-D centroid coordinate data set
  centroid_coord_high_D <- centroid_coord_high_D |>
    dplyr::mutate(distance = distance_vec) |>
    dplyr::mutate(hb_id = avg_df_training$hb_id)

  ## Sort by distance and obtain the centroid which is nearest
  predict_centroid_coord_high_D <- centroid_coord_high_D |>
    dplyr::arrange(distance) |>
    dplyr::filter(dplyr::row_number() == 1)

  ## Rename columns
  names(predict_centroid_coord_high_D)[1:(NCOL(test_data_point) - 1)] <- paste0("C_", names(predict_centroid_coord_high_D)[1:(NCOL(test_data_point) - 1)])

  ## Obtain 2D coordinate of the nearest high-D centroid
  predict_centroid_coord_2D <- centroid_df_training |>
    dplyr::filter(hexID %in% predict_centroid_coord_high_D$hb_id) |>
    dplyr::select(x, y)

  ## Rename columns
  names(predict_centroid_coord_2D) <- paste0("C_", names(predict_centroid_coord_2D))

  ## Combine high-D and 2D coordinate
  predict_centroid_coord_all <- dplyr::bind_cols(predict_centroid_coord_high_D, predict_centroid_coord_2D)

  ## Combine all
  predict_coord_test <- dplyr::bind_rows(predict_coord_test, predict_centroid_coord_all)


}

### Measure

UMAP_data_predict <- predict_coord_test |>
  dplyr::select("C_x", "C_y") |>
  dplyr::mutate(ID = test_data$ID)

names(UMAP_data_predict)[1:2] <- c("UMAP1", "UMAP2")

pred_df_test_object <- predict_hex_id(training_data = training_data, nldr_df = UMAP_data, nldr_df_test = UMAP_data_predict, num_bins = num_bins_curve, shape_val = shape_value_curve, x = "UMAP1", y = "UMAP2", col_start = "x")
pred_df_test <- pred_df_training_object$pred_data
centroid_df_test <- pred_df_training_object$df_bin_centroids
avg_df_test <- pred_df_training_object$df_bin

eval_df_test <- generate_eval_df(data = data, prediction_df = pred_df_test, df_bin_centroids = centroid_df_test, df_bin = avg_df_test, num_bins = num_bins_curve, col_start = "x")
