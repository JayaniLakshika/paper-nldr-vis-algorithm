library(readr)
library(dplyr)


set.seed(20240110)

source("quollr_code.R", local = TRUE)

data <- read_rds("data/zeisel/zeisel.rds")
training_data <- read_rds("data/zeisel/zeisel_training.rds")
test_data <- read_rds("data/zeisel/zeisel_test.rds")

tSNE_data <- read_rds("data/zeisel/zeisel_tsne_30.rds")
UMAP_data <- read_rds("data/zeisel/zeisel_umap.rds")
PHATE_data <- read_rds("data/zeisel/zeisel_phate.rds")
TriMAP_data <- read_rds("data/zeisel/zeisel_trimap.rds")
PaCMAP_data <- read_rds("data/zeisel/zeisel_pacmap.rds")


shape_value_curve <- calculate_effective_shape_value(.data = UMAP_data,
                                                     x = UMAP1, y = UMAP2)

num_bins_x <- 6 ## Number of bins along the x-axis

benchmark_rm_hex_vec <- sample(seq(0, 0.3, by=0.01), 10) |> sort()

vec <- stats::setNames(rep("", 7), c("number_of_bins", "number_of_observations", "total_error", "totol_error_method_2", "totol_error_method_3", "total_mse", "benchmark_rm_hex"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(benchmark_rm_hex_vec)) {

  pred_df_training_object <- predict_hex_id(training_data = training_data, nldr_df = UMAP_data, nldr_df_test = UMAP_data, num_bins = num_bins_x, shape_val = shape_value_curve, x = "UMAP1", y = "UMAP2", col_start = "PC")
  pred_df_training <- pred_df_training_object$pred_data
  centroid_df_training_all <- pred_df_training_object$df_bin_centroids

  hexbin_data_object_training <- pred_df_training_object$hexbin_data_object

  ## Identify bins with low-density
  identify_rm_bins <- find_low_density_hexagons(centroid_df_training_all, num_bins_x, benchmark_rm_hex = benchmark_rm_hex_vec[i])

  centroid_df_training <- centroid_df_training_all |>
    filter(!(hexID %in% identify_rm_bins))

  ## Add hexbin Id to 2D embeddings
  UMAP_data_with_hb_id <- UMAP_data |>
    mutate(hb_id = hexbin_data_object_training$hb_data@cID)

  ## To generate a data set with high-D and 2D training data
  df_all <- dplyr::bind_cols(training_data |> dplyr::select(-ID), UMAP_data_with_hb_id)

  ## Averaged on high-D
  avg_df_training <- avg_highD_data(.data = df_all, column_start_text = "PC")

  eval_df_training <- generate_eval_df(data = data, prediction_df = pred_df_training,
                                       df_bin_centroids_all = centroid_df_training_all, df_bin = avg_df_training,
                                       num_bins = num_bins_x, col_start = "PC", df_bin_centroids = centroid_df_training)

  eval_df_training <- eval_df_training |>
    dplyr::mutate(benchmark_rm_hex = benchmark_rm_hex_vec[i])

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)


}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

MSE_df_2 <- eval_data_training |>
  dplyr::mutate(method = "UMAP")

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
