library(readr)
library(dplyr)
library(ggplot2)

set.seed(20240110)

source("quollr_code.R", local = TRUE)

data <- read_rds("data/s_curve/s_curve.rds")
training_data <- read_rds("data/s_curve/s_curve_training.rds")
test_data <- read_rds("data/s_curve/s_curve_test.rds")

UMAP_data <- read_rds("data/s_curve/s_curve_umap.rds")

shape_value_curve <- calculate_effective_shape_value(.data = UMAP_data,
                                                     x = UMAP1, y = UMAP2)

num_bins_x <- 8 ## Number of bins along the x-axis (looking at mse plot)

vec <- stats::setNames(rep("", 4), c("number_of_bins", "number_of_observations",
                                     "total_error", "total_mse"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

model_object <- fit_high_d_model(training_data = training_data,
                                 nldr_df_with_id = UMAP_data, x = "UMAP1", y = "UMAP2",
                                 num_bins_x = num_bins_x, shape_val = shape_value_curve,
                                 is_bin_centroid = TRUE,
                                 is_rm_lwd_hex = FALSE,
                                 benchmark_to_rm_lwd_hex = NA,
                                 is_avg_high_d = TRUE, column_start_text = "x")

centroid_df_training_all <- model_object$df_bin_centroids
avg_df_training <- model_object$df_bin

benchmark_rm_hex_vec <- seq(0, 1, by=0.1)
benchmark_rm_hex_vec <- append(benchmark_rm_hex_vec[1:10],0.99)


for (i in 1:length(benchmark_rm_hex_vec)) {

  centroid_df_training <- centroid_df_training_all |>
    filter(std_counts > benchmark_rm_hex_vec[i])

  pred_df_training <- predict_2d_embeddings(test_data = training_data,
                                            df_bin_centroids = centroid_df_training,
                                            df_bin = avg_df_training, type_NLDR = "UMAP")

  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training,
                                       df_bin_centroids = centroid_df_training,
                                       df_bin = avg_df_training, col_start = "x")

  eval_df_training <- eval_df_training |>
    dplyr::mutate(benchmark_rm_hex = benchmark_rm_hex_vec[i])

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

MSE_df_2 <- eval_data_training |>
  dplyr::mutate(method = "UMAP")

write_rds(MSE_df_2, "data/s_curve/s_curve_summary_lwd_umap.rds")

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
