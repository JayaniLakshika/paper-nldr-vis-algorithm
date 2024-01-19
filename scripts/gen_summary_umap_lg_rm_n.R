library(readr)
library(dplyr)

set.seed(20240110)

source("quollr_code.R", local = TRUE)
source("nldr_code.R", local = TRUE)

## Import data
df_s_curve <- read_rds("data/s_curve/s_curve.rds")
training_data_s_curve <- read_rds("data/s_curve/s_curve_training.rds")

UMAP_s_curve <- read_rds("data/s_curve/s_curve_umap.rds")

num_bins_s_curve <- 6

shape_value_s_curve <- calculate_effective_shape_value(.data = UMAP_s_curve,
                                                      x = UMAP1, y = UMAP2) ## 1.259938

## To extract bin centroids
hexbin_data_object_s_curve <- extract_hexbin_centroids(nldr_df = UMAP_s_curve,
                                                      num_bins = num_bins_s_curve,
                                                      shape_val = shape_value_s_curve, x = UMAP1, y = UMAP2)

df_bin_centroids_s_curve <- hexbin_data_object_s_curve$hexdf_data

UMAP_s_curve_with_hb_id <- UMAP_s_curve |>
  dplyr::mutate(hb_id = hexbin_data_object_s_curve$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_s_curve <- dplyr::bind_cols(training_data_s_curve |> dplyr::select(-ID), UMAP_s_curve_with_hb_id)

## Averaged on high-D
df_bin_s_curve <- avg_highD_data(.data = df_all_s_curve, column_start_text = "x")

## Triangulate bin centroids
tr1_object_s_curve <- triangulate_bin_centroids(df_bin_centroids_s_curve, x, y)
tr_from_to_df_s_curve <- generate_edge_info(triangular_object = tr1_object_s_curve)

## Compute 2D distances
distance_s_curve <- cal_2D_dist(.data = tr_from_to_df_s_curve)

## To find the benchmark value
benchmark <- find_benchmark_value(.data = distance_s_curve, distance_col = distance)

benchmark_dist_vec <- distance_s_curve$distance |> round(3) |> unique() |> sort()
#benchmark_dist_vec <- seq(min(distance_s_curve$distance |> round(3) |> unique()), max(distance_s_curve$distance |> round(3) |> unique()), 1)

vec <- stats::setNames(rep("", 2), c("benchmark_rm_lg", "total_error"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for(i in 1:length(benchmark_dist_vec)) {

  pred_df_training_object <- predict_hex_id(training_data = training_data_s_curve, nldr_df = UMAP_s_curve, nldr_df_test = UMAP_s_curve, num_bins = num_bins_s_curve, shape_val = shape_value_s_curve, x = "UMAP1", y = "UMAP2", col_start = "x")
  pred_df_training <- pred_df_training_object$pred_data
  centroid_df_training <- pred_df_training_object$df_bin_centroids
  avg_df_training <- pred_df_training_object$df_bin

  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, num_bins = num_bins_vec[i], col_start = "x")

  eval_df_training <- eval_df_training |>
    mutate(benchmark_rm_lg = benchmark_dist_vec[i])

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}

ggplot(eval_data_training, aes(x = benchmark_rm_lg,
                                                               y = total_error
)) +
  geom_point() +
  geom_line()

eval_data_training <- eval_data_training |>
  dplyr::mutate(method = "UMAP")

write_rds(eval_data_training, "data/s_curve/s_curve_summary_lg_umap.rds")
