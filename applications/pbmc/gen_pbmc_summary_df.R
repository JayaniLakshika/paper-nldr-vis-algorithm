library(readr)
library(dplyr)

set.seed(20240110)

source("quollr_code.R", local = TRUE)

data <- read_csv("pbmc/pbmc.rds")
training_data <- read_csv("pbmc/pbmc_training.rds")
test_data <- read_csv("pbmc/pbmc_test.rds")

tSNE_data <- read_csv("pbmc/pbmc_tsne_30.rds")
UMAP_data <- read_csv("pbmc/pbmc_umap.rds")
PHATE_data <- read_csv("pbmc/pbmc_phate.rds")
TriMAP_data <- read_csv("pbmc/pbmc_trimap.rds")
PaCMAP_data <- read_csv("pbmc/pbmc_pacmap.rds")


## tSNE

shape_value_curve <- calculate_effective_shape_value(.data = tSNE_data,
                                                     x = tSNE1, y = tSNE2)

num_bins_vec <- 1:10 ## Number of bins along the x-axis

vec <- stats::setNames(rep("", 6), c("number_of_bins", "number_of_observations", "total_error", "totol_error_method_2", "totol_error_method_3", "total_mse"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(num_bins_vec)) {

  pred_df_training_object <- predict_hex_id(training_data = training_data, nldr_df = tSNE_data, nldr_df_test = tSNE_data, num_bins = num_bins_vec[i], shape_val = shape_value_curve, x = "tSNE1", y = "tSNE2", col_start = "x")
  pred_df_training <- pred_df_training_object$pred_data
  centroid_df_training <- pred_df_training_object$df_bin_centroids
  avg_df_training <- pred_df_training_object$df_bin

  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, num_bins = num_bins_vec[i], col_start = "x")

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

MSE_df_1 <- eval_data_training |>
  dplyr::mutate(method = "tSNE")


shape_value_curve <- calculate_effective_shape_value(.data = UMAP_data,
                                                     x = UMAP1, y = UMAP2)

num_bins_vec <- 1:10 ## Number of bins along the x-axis

vec <- stats::setNames(rep("", 6), c("number_of_bins", "number_of_observations", "total_error", "totol_error_method_2", "totol_error_method_3", "total_mse"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(num_bins_vec)) {

  pred_df_training_object <- predict_hex_id(training_data = training_data, nldr_df = UMAP_data, nldr_df_test = UMAP_data, num_bins = num_bins_vec[i], shape_val = shape_value_curve, x = "UMAP1", y = "UMAP2", col_start = "x")
  pred_df_training <- pred_df_training_object$pred_data
  centroid_df_training <- pred_df_training_object$df_bin_centroids
  avg_df_training <- pred_df_training_object$df_bin

  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, num_bins = num_bins_vec[i], col_start = "x")

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

MSE_df_2 <- eval_data_training |>
  dplyr::mutate(method = "UMAP")



shape_value_curve <- calculate_effective_shape_value(.data = PHATE_data,
                                                     x = PHATE1, y = PHATE2)

num_bins_vec <- 1:10 ## Number of bins along the x-axis

vec <- stats::setNames(rep("", 6), c("number_of_bins", "number_of_observations", "total_error", "totol_error_method_2", "totol_error_method_3", "total_mse"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(num_bins_vec)) {

  pred_df_training_object <- predict_hex_id(training_data = training_data, nldr_df = PHATE_data, nldr_df_test = PHATE_data, num_bins = num_bins_vec[i], shape_val = shape_value_curve, x = "PHATE1", y = "PHATE2", col_start = "x")
  pred_df_training <- pred_df_training_object$pred_data
  centroid_df_training <- pred_df_training_object$df_bin_centroids
  avg_df_training <- pred_df_training_object$df_bin

  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, num_bins = num_bins_vec[i], col_start = "x")

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

MSE_df_3 <- eval_data_training |>
  dplyr::mutate(method = "PHATE")


## TriMAP

## Prediction

shape_value_curve <- calculate_effective_shape_value(.data = TriMAP_data,
                                                     x = TriMAP1, y = TriMAP2)

num_bins_vec <- 1:10 ## Number of bins along the x-axis

vec <- stats::setNames(rep("", 6), c("number_of_bins", "number_of_observations", "total_error", "totol_error_method_2", "totol_error_method_3", "total_mse"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(num_bins_vec)) {

  pred_df_training_object <- predict_hex_id(training_data = training_data, nldr_df = TriMAP_data, nldr_df_test = TriMAP_data, num_bins = num_bins_vec[i], shape_val = shape_value_curve, x = "TriMAP1", y = "TriMAP2", col_start = "x")
  pred_df_training <- pred_df_training_object$pred_data
  centroid_df_training <- pred_df_training_object$df_bin_centroids
  avg_df_training <- pred_df_training_object$df_bin

  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, num_bins = num_bins_vec[i], col_start = "x")

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

MSE_df_4 <- eval_data_training |>
  dplyr::mutate(method = "TriMAP")



shape_value_curve <- calculate_effective_shape_value(.data = PaCMAP_data,
                                                     x = PaCMAP1, y = PaCMAP2)

num_bins_vec <- 1:10 ## Number of bins along the x-axis

vec <- stats::setNames(rep("", 6), c("number_of_bins", "number_of_observations", "total_error", "totol_error_method_2", "totol_error_method_3", "total_mse"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(num_bins_vec)) {

  pred_df_training_object <- predict_hex_id(training_data = training_data, nldr_df = PaCMAP_data, nldr_df_test = PaCMAP_data, num_bins = num_bins_vec[i], shape_val = shape_value_curve, x = "PaCMAP1", y = "PaCMAP2", col_start = "x")
  pred_df_training <- pred_df_training_object$pred_data
  centroid_df_training <- pred_df_training_object$df_bin_centroids
  avg_df_training <- pred_df_training_object$df_bin

  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, num_bins = num_bins_vec[i], col_start = "x")

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

MSE_df_5 <- eval_data_training |>
  dplyr::mutate(method = "PaCMAP")


MSE_df <- dplyr::bind_rows(MSE_df_1, MSE_df_2, MSE_df_3, MSE_df_4, MSE_df_5)
