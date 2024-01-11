library(readr)
library(dplyr)
library(ggplot2)
library(colorspace)

set.seed(20240110)

source("quollr_code.R", local = TRUE)
source("nldr_code.R", local = TRUE)

data <- read_rds("data/s_curve/s_curve.rds")
training_data <- read_rds("data/s_curve/s_curve_training.rds")
UMAP_data <- read_rds("data/s_curve/s_curve_umap_15.rds")

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


## To draw with AIC
ggplot(MSE_df_2 |> dplyr::filter(data_type == "training"), aes(x = number_of_bins,
                                                               y = total_error,
                                                               color = method
)) +
  geom_point() +
  geom_line() +
  #geom_vline(xintercept = NROW(full_grid_with_hexbin_id)) +
  #annotate("text", x= (NROW(full_grid_with_hexbin_id) - 10), y=-5000, label=paste0("effective number of bins = ", as.character(NROW(full_grid_with_hexbin_id))), angle=90) +
  theme_light() +
  theme(legend.title = element_blank(), plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7)) +
  scale_color_discrete_qualitative() +
  ylab("AIC") +
  xlab("Total number of bins")
## Effective number of bins along x-axis

ggplot(MSE_df_2, aes(x = number_of_bins,
                     y = total_mse,
                     color = method
)) +
  geom_point() +
  geom_line() +
  theme_light() +
  theme(legend.title = element_blank(), plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7)) +
  # geom_vline(xintercept = NROW(full_grid_with_hexbin_id)) +
  # annotate("text", x= (NROW(full_grid_with_hexbin_id) - 10), y=0.25, label=paste0("effective number of bins = ", as.character(NROW(full_grid_with_hexbin_id))), angle=90) +
  scale_color_discrete_qualitative() +
  ylab("MSE") +
  xlab("Total number of bins")
