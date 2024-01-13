library(readr)
library(dplyr)
library(ggplot2)

set.seed(20240110)

source("quollr_code.R", local = TRUE)

data <- read_rds("data/s_curve/s_curve.rds")
training_data <- read_rds("data/s_curve/s_curve_training.rds")
test_data <- read_rds("data/s_curve/s_curve_test.rds")

PHATE_data <- read_rds("data/s_curve/s_curve_phate.rds")

shape_value_curve <- calculate_effective_shape_value(.data = PHATE_data,
                                                     x = PHATE1, y = PHATE2)

num_bins_x <- 5 ## Number of bins along the x-axis (looking at mse plot)

benchmark_rm_hex_vec <- sample(seq(0, 0.5, by=0.001), 20) |> sort()

vec <- stats::setNames(rep("", 4), c("number_of_bins", "number_of_observations", "total_error", "total_mse"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

pred_df_training_object <- predict_hex_id(training_data = training_data, nldr_df = PHATE_data, nldr_df_test = PHATE_data, num_bins = num_bins_x, shape_val = shape_value_curve, x = "PHATE1", y = "PHATE2", col_start = "x")
pred_df_training <- pred_df_training_object$pred_data
centroid_df_training_all <- pred_df_training_object$df_bin_centroids

hexbin_data_object_training <- pred_df_training_object$hexbin_data_object

for (i in 1:length(benchmark_rm_hex_vec)) {

  ## Identify bins with low-density
  identify_rm_bins <- find_low_density_hexagons(centroid_df_training_all, num_bins_x, benchmark_rm_hex = benchmark_rm_hex_vec[i])

  centroid_df_training <- centroid_df_training_all |>
    filter(!(hexID %in% identify_rm_bins))

  ## Add hexbin Id to 2D embeddings
  PHATE_data_with_hb_id <- PHATE_data |>
    mutate(hb_id = hexbin_data_object_training$hb_data@cID)

  ## To generate a data set with high-D and 2D training data
  df_all <- dplyr::bind_cols(training_data |> dplyr::select(-ID), PHATE_data_with_hb_id)

  ## Averaged on high-D
  avg_df_training <- avg_highD_data(.data = df_all, column_start_text = "x")

  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training,
                                       df_bin_centroids_all = centroid_df_training_all, df_bin = avg_df_training,
                                       num_bins = num_bins_x, col_start = "x", df_bin_centroids = centroid_df_training, rm_lwd_hex = TRUE)

  eval_df_training <- eval_df_training |>
    dplyr::mutate(benchmark_rm_hex = benchmark_rm_hex_vec[i])

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

MSE_df_2 <- eval_data_training |>
  dplyr::mutate(method = "PHATE")

write_rds(MSE_df_2, "data/s_curve/s_curve_summary_lwd_phate.rds")

## To draw with AIC
ggplot(MSE_df_2 |> dplyr::filter(data_type == "training"), aes(x = benchmark_rm_hex,
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
  scale_colour_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  ylab("AIC") +
  xlab("Benchmark value to remove low-density hexagons")
## Effective number of bins along x-axis

ggplot(MSE_df_2, aes(x = benchmark_rm_hex,
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
  scale_colour_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  ylab("MSE") +
  xlab("Benchmark value to remove low-density hexagons")
