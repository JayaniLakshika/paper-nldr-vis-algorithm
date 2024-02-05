library(readr)

set.seed(20240110)

source("quollr_code.R", local = TRUE)

## UMAP
## Prediction
UMAP_data <- read_rds(file = "data/s_curve/s_curve_umap.rds")
hex_full_count_df <- read_rds("data/s_curve/s_curve_hex_8.rds")

shape_val_umap_s_curve <- calculate_effective_shape_value(.data = UMAP_data,
                                                          x = UMAP1, y = UMAP2)

hexbin_data_object <- extract_hexbin_centroids(UMAP_data, 8, shape_val_umap_s_curve)

## Add hexbin Id to 2D embeddings
UMAP_data_with_hb_id <- UMAP_data |>
  dplyr::mutate(hb_id = hexbin_data_object$hb_data@cID)

shift_amount_vec <- append(seq(0,0.537285, 0.1) * (-1), seq(0,0.537285, 0.1)) |> unique() ## Shift amount

vec <- stats::setNames(rep("", 6), c("number_of_bins", "number_of_observations", "total_error", "total_mse", "num_bins_x", "shift"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(shift_amount_vec)) { #length(shift_amount_vec)

  shifted_object <- extract_coord_of_shifted_hex_grid(nldr_data_with_hb_id = UMAP_data_with_hb_id, num_bins_x = 8, hex_full_count_df, shift_x = shift_amount_vec[i], shift_y = shift_amount_vec[i])

  centroid_df_training <- shifted_object$hex_full_count_df_new |>
    dplyr::filter(!is.na(counts)) |>
    dplyr::select(c_x, c_y, hexID, counts, std_counts) |>
    dplyr::rename(c("x" = "c_x",
                    "y" = "c_y")) |>
    dplyr::distinct()

  # pred_df_training <- predict_2d_embeddings(test_data = training_data, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, type_NLDR = "tSNE")

  ## Add hexbin Id to 2D embeddings
  UMAP_data_with_hb_id_n <- shifted_object$nldr_df_with_new_hexID

  ## To generate a data set with high-D and 2D training data
  df_all_umap_s_curve <- dplyr::bind_cols(training_data |> dplyr::select(-ID), UMAP_data_with_hb_id_n)

  ## Averaged on high-D
  avg_df_training <- avg_highD_data(.data = df_all_umap_s_curve)

  pred_df_training <- predict_hex_id(df_bin_centroids = centroid_df_training, nldr_df_test = UMAP_data, x = "UMAP1", y = "UMAP2")

  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, num_bins = 8, col_start = "x")

  eval_df_training <- eval_df_training |>
    mutate(num_bins_x = 8,
           shift = shift_amount_vec[i])

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

MSE_df_2 <- eval_data_training |>
  dplyr::mutate(method = "UMAP")

write_rds(MSE_df_2, "data/s_curve/s_curve_umap_shift_summary.rds")
