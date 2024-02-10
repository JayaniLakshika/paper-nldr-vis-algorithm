```{r}
#| warning: false
#| echo: false



## tSNE

shape_value_curve <- calculate_effective_shape_value(.data = tSNE_data,
                                                     x = tSNE1, y = tSNE2)

num_bins_vec <- 1:18 ## Number of bins along the x-axis

vec <- stats::setNames(rep("", 5), c("number_of_bins", "number_of_observations", "total_error", "total_mse", "num_bins_x"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(num_bins_vec)) {

  model_object <- fit_high_d_model(training_data = training_data, nldr_df_with_id = tSNE_data, x = "tSNE1", y = "tSNE2", num_bins_x = num_bins_vec[i], shape_val = shape_value_curve,
                                   is_bin_centroid = TRUE,
                                   is_rm_lwd_hex = FALSE,
                                   benchmark_to_rm_lwd_hex = NA,
                                   is_avg_high_d = TRUE, column_start_text = "x")

  centroid_df_training <- model_object$df_bin_centroids
  avg_df_training <- model_object$df_bin

  # pred_df_training <- predict_2d_embeddings(test_data = training_data, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, type_NLDR = "tSNE")

  pred_df_training <- predict_hex_id(df_bin_centroids = centroid_df_training, nldr_df_test = tSNE_data, x = "tSNE1", y = "tSNE2")


  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, col_start = "x")

  eval_df_training <- eval_df_training |>
    mutate(num_bins_x = num_bins_vec[i])

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

MSE_df_1 <- eval_data_training |>
  dplyr::mutate(method = "tSNE")
```
```{r}
#| warning: false
#| echo: false
## PAHTE
## Prediction

shape_value_curve <- calculate_effective_shape_value(.data = PHATE_data,
                                                     x = PHATE1, y = PHATE2)

num_bins_vec <- 1:18 ## Number of bins along the x-axis

vec <- stats::setNames(rep("", 5), c("number_of_bins", "number_of_observations", "total_error", "total_mse", "num_bins_x"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(num_bins_vec)) {

  model_object <- fit_high_d_model(training_data = training_data, nldr_df_with_id = PHATE_data, x = "PHATE1", y = "PHATE2", num_bins_x = num_bins_vec[i], shape_val = shape_value_curve,
                                   is_bin_centroid = TRUE,
                                   is_rm_lwd_hex = FALSE,
                                   benchmark_to_rm_lwd_hex = NA,
                                   is_avg_high_d = TRUE, column_start_text = "x")

  centroid_df_training <- model_object$df_bin_centroids
  avg_df_training <- model_object$df_bin

  # pred_df_training <- predict_2d_embeddings(test_data = training_data, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, type_NLDR = "tSNE")

  pred_df_training <- predict_hex_id(df_bin_centroids = centroid_df_training, nldr_df_test = PHATE_data, x = "PHATE1", y = "PHATE2")


  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, col_start = "x")

  eval_df_training <- eval_df_training |>
    mutate(num_bins_x = num_bins_vec[i])

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

MSE_df_3 <- eval_data_training |>
  dplyr::mutate(method = "PHATE")

```

```{r}
#| warning: false
#| echo: false

## TriMAP

## Prediction

shape_value_curve <- calculate_effective_shape_value(.data = TriMAP_data,
                                                     x = TriMAP1, y = TriMAP2)

num_bins_vec <- 1:10 ## Number of bins along the x-axis

vec <- stats::setNames(rep("", 5), c("number_of_bins", "number_of_observations", "total_error", "total_mse", "num_bins_x"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(num_bins_vec)) {

  model_object <- fit_high_d_model(training_data = training_data, nldr_df_with_id = TriMAP_data, x = "TriMAP1", y = "TriMAP2", num_bins_x = num_bins_vec[i], shape_val = shape_value_curve,
                                   is_bin_centroid = TRUE,
                                   is_rm_lwd_hex = FALSE,
                                   benchmark_to_rm_lwd_hex = NA,
                                   is_avg_high_d = TRUE, column_start_text = "x")

  centroid_df_training <- model_object$df_bin_centroids
  avg_df_training <- model_object$df_bin

  # pred_df_training <- predict_2d_embeddings(test_data = training_data, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, type_NLDR = "tSNE")

  pred_df_training <- predict_hex_id(df_bin_centroids = centroid_df_training, nldr_df_test = TriMAP_data, x = "TriMAP1", y = "TriMAP2")


  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, col_start = "x")

  eval_df_training <- eval_df_training |>
    mutate(num_bins_x = num_bins_vec[i])

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

MSE_df_4 <- eval_data_training |>
  dplyr::mutate(method = "TriMAP")

```

```{r}
#| warning: false
#| echo: false

## PaCMAP

## Prediction

shape_value_curve <- calculate_effective_shape_value(.data = PaCMAP_data,
                                                     x = PaCMAP1, y = PaCMAP2)

num_bins_vec <- 1:18 ## Number of bins along the x-axis

vec <- stats::setNames(rep("", 5), c("number_of_bins", "number_of_observations", "total_error", "total_mse", "num_bins_x"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(num_bins_vec)) {

  model_object <- fit_high_d_model(training_data = training_data, nldr_df_with_id = PaCMAP_data, x = "PaCMAP1", y = "PaCMAP2", num_bins_x = num_bins_vec[i], shape_val = shape_value_curve,
                                   is_bin_centroid = TRUE,
                                   is_rm_lwd_hex = FALSE,
                                   benchmark_to_rm_lwd_hex = NA,
                                   is_avg_high_d = TRUE, column_start_text = "x")

  centroid_df_training <- model_object$df_bin_centroids
  avg_df_training <- model_object$df_bin

  # pred_df_training <- predict_2d_embeddings(test_data = training_data, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, type_NLDR = "tSNE")

  pred_df_training <- predict_hex_id(df_bin_centroids = centroid_df_training, nldr_df_test = PaCMAP_data, x = "PaCMAP1", y = "PaCMAP2")


  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, col_start = "x")

  eval_df_training <- eval_df_training |>
    mutate(num_bins_x = num_bins_vec[i])

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

MSE_df_5 <- eval_data_training |>
  dplyr::mutate(method = "PaCMAP")

```
