####### S-curve example with 5-Fold cross validation #####################
library(readr)

set.seed(20240110)

source("quollr_code.R", local = TRUE)

## Import data
data <- read_rds("data/s_curve/s_curve.rds")
training_data <- read_rds("data/s_curve/s_curve_training.rds")
UMAP_s_curve <- read_rds("data/s_curve/s_curve_umap.rds")

shape_val_umap_s_curve <- calculate_effective_shape_value(.data = UMAP_s_curve,
                                                          x = UMAP1, y = UMAP2)

# Create indices for 5-fold cross-validation
num_obs <- NROW(UMAP_s_curve)
fold_size <- floor(num_obs / 5)

# Initialize variables to store performance metrics
mse_list <- numeric(5)

num_bins_vec <- 1:13

vec <- stats::setNames(rep("", 5), c("number_of_bins", "number_of_observations",
                                     "total_error", "total_mse", "num_bins_x"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

# Perform 5-fold cross-validation
for (fold in 1:5) {
  # Define the indices for the current fold
  start_index <- (fold - 1) * fold_size + 1
  end_index <- fold * fold_size

  # Split the data into training and validation sets
  X_train <- UMAP_s_curve[-(start_index:end_index), ]
  X_val <- UMAP_s_curve[start_index:end_index, ]

  for (i in 1:length(num_bins_vec)) {
    # Fit the model for training
    model_object <- fit_high_d_model(training_data = training_data, nldr_df_with_id = UMAP_s_curve,
                                     x = "UMAP1", y = "UMAP2", num_bins_x = num_bins_vec[i], shape_val = shape_val_umap_s_curve,
                                     is_bin_centroid = TRUE,
                                     is_rm_lwd_hex = FALSE,
                                     benchmark_to_rm_lwd_hex = NA,
                                     is_avg_high_d = TRUE, column_start_text = "x")

    centroid_df_training <- model_object$df_bin_centroids
    avg_df_training <- model_object$df_bin

    # Make predictions on the validation set
    pred_df_training <- predict_hex_id(df_bin_centroids = centroid_df_training,
                                       nldr_df_test = X_val, x = "UMAP1", y = "UMAP2")

    # Calculate mean squared error (MSE) for the fold
    eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training,
                                         df_bin_centroids = centroid_df_training, df_bin = avg_df_training,
                                         num_bins = num_bins_vec[i], col_start = "x")

    eval_df_training <- eval_df_training |>
      mutate(num_bins_x = num_bins_vec[i],
             fold = fold)

    eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)

  }



}

eval_data_training_n <- eval_data_training |>
  group_by(number_of_bins, num_bins_x, number_of_observations) |>
  summarise(mean_mse = mean(total_mse),
            mean_aic = mean(total_error), .groups = "drop")

## To draw with AIC
aic_plot <- ggplot(eval_data_training_n, aes(x = number_of_bins, y = mean_aic
)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 144, linetype="dashed",
             color = "red", size=0.5) +
  #geom_vline(xintercept = NROW(full_grid_with_hexbin_id)) +
  #annotate("text", x= (NROW(full_grid_with_hexbin_id) - 10), y=-5000, label=paste0("effective number of bins = ", as.character(NROW(full_grid_with_hexbin_id))), angle=90) +
  theme_light() +
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7)) +
  scale_colour_manual(values = c("#377eb8", "#e41a1c", "#4daf4a", "#984ea3", "#ff7f00")) +
  ylab("AIC") +
  xlab("total number of bins")
## Effective number of bins along x-axis

mse_plot <- ggplot(eval_data_training_n, aes(x = number_of_bins,
                                                                  y = mean_mse
)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 144, linetype="dashed",
             color = "red", size=0.5) +
  theme_light() +
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7)) +
  # geom_vline(xintercept = NROW(full_grid_with_hexbin_id)) +
  # annotate("text", x= (NROW(full_grid_with_hexbin_id) - 10), y=0.25, label=paste0("effective number of bins = ", as.character(NROW(full_grid_with_hexbin_id))), angle=90) +
  scale_colour_manual(values = c("#377eb8", "#e41a1c", "#4daf4a", "#984ea3", "#ff7f00")) +
  ylab("MSE") +
  xlab("total number of bins")




