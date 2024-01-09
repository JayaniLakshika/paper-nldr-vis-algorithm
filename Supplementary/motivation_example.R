library(dplyr)
library(snedata)
library(ggflowchart)
library(purrr) ## map function
library(gridExtra) ## for grid.arrange
library(rsample)
library(DT)
library(ggbeeswarm)
library(ggplot2)
library(readr)

library(Rtsne)
library(umap)
library(phateR)
library(reticulate)
library(patchwork)

library(grid)
set.seed(20230531)

source(paste0(here::here(), "/paper-nldr-vis-algorithm/quollr_code.R"))
source(paste0(here::here(), "/paper-nldr-vis-algorithm/nldr_code.R"))

sample_size <- 5000
cluster_size <- sample_size/5
df1 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.05), y=rnorm(cluster_size, mean = 0, sd = 0.05), z=rnorm(cluster_size, mean = 0, sd = 0.05), w=rnorm(cluster_size, mean = 0, sd = 0.05))

df2 <- tibble::tibble(x=rnorm(cluster_size, mean = 1, sd = 0.05), y=rnorm(cluster_size, mean = 0, sd = 0.05), z=rnorm(cluster_size, mean = 0, sd = 0.05), w=rnorm(cluster_size, mean = 0, sd = 0.05))

df3 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.05), y=rnorm(cluster_size, mean = 1, sd = 0.05), z=rnorm(cluster_size, mean = 0, sd = 0.05), w=rnorm(cluster_size, mean = 0, sd = 0.05))

df4 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.05), y=rnorm(cluster_size, mean = 0, sd = 0.05), z=rnorm(cluster_size, mean = 1, sd = 0.05), w=rnorm(cluster_size, mean = 0, sd = 0.05))

df5 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.05), y=rnorm(cluster_size, mean = 0, sd = 0.05), z=rnorm(cluster_size, mean = 0, sd = 0.05), w=rnorm(cluster_size, mean = 1, sd = 0.05))

df_2 <- bind_rows(df1, df2, df3, df4, df5)
df_2 <- df_2 %>%
  rename(x1 = x, x2 = y, x3 = z, x4 = w)

df_2 <- df_2 |>
  mutate(ID = row_number())

data_split_sp <- initial_split(df_2)
training_data_5 <- training(data_split_sp) |>
  arrange(ID)
test_data_5 <- testing(data_split_sp) |>
  arrange(ID)

### tSNE
tSNE_data_gau <- Fit_tSNE(training_data_5 |> dplyr::select(-ID), opt_perplexity = calculate_effective_perplexity(training_data_5), with_seed = 20230531)

plot_gau <- plot_tSNE_2D(tSNE_data_gau) + #ggtitle("(a)") +
  theme_linedraw() +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())



cell_area <- 1

num_bins_tsne_gau <- calculate_effective_x_bins(.data = tSNE_data_gau, x = tSNE1,
                                                cell_area = 1)
shape_val_tsne_gau <- calculate_effective_shape_value(.data = tSNE_data_gau,
                                                      x = tSNE1, y = tSNE2)

## To extract bin centroids
hexbin_data_object_tsne_gau <- extract_hexbin_centroids(nldr_df = tSNE_data_gau, num_bins = num_bins_tsne_gau, shape_val = shape_val_tsne_gau, x = tSNE1, y = tSNE2)

df_bin_centroids_tsne_gau <- hexbin_data_object_tsne_gau$hexdf_data

# min_std_cell_threshold3 <- 0.01
#
# df_bin_centroids_pbmc <- df_bin_centroids_pbmc |>
#   dplyr::mutate(stand_cell_count = counts/max(counts)) |>
#   dplyr::filter(stand_cell_count > min_std_cell_threshold3)

tSNE_data_with_hb_id_gau <- tSNE_data_gau |>
  dplyr::mutate(hb_id = hexbin_data_object_tsne_gau$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_tsne_gau <- dplyr::bind_cols(training_data_5 |> dplyr::select(-ID), tSNE_data_with_hb_id_gau)

## Averaged on high-D
df_bin_tsne_gau <- avg_highD_data(.data = df_all_tsne_gau)

## Triangulate bin centroids
tr1_object_tsne_gau <- triangulate_bin_centroids(df_bin_centroids_tsne_gau, x, y)
tr_from_to_df_tsne_gau <- generate_edge_info(triangular_object = tr1_object_tsne_gau)

## Compute 2D distances
distance_tsne_gau <- cal_2D_dist(.data = tr_from_to_df_tsne_gau)

## To find the benchmark value
benchmark_tsne_gau <- find_benchmark_value(.data = distance_tsne_gau, distance_col = distance)
#benchmark_tsne <- 5

tour_tsne_gau <- show_langevitour(df_all_tsne_gau, df_bin_tsne_gau, df_bin_centroids_tsne_gau, benchmark_value = benchmark_tsne_gau, distance = distance_tsne_gau, distance_col = distance)

## Prediction

shape_value <- calculate_effective_shape_value(.data = tSNE_data_gau,
                                               x = tSNE1, y = tSNE2)

num_bins_vec <- 1:10 ## Number of bins along the x-axis

vec <- stats::setNames(rep("", 6), c("number_of_bins", "number_of_observations", "total_error", "totol_error_method_2", "totol_error_method_3", "total_mse"))  ## Define column names

eval_data_test <- dplyr::bind_rows(vec)[0, ]
eval_data_test <- eval_data_test |>
  dplyr::mutate_if(is.character, as.numeric)

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(num_bins_vec)) {

  pred_df_training_object <- predict_hex_id(training_data = training_data_5, nldr_df = tSNE_data_gau, nldr_df_test = tSNE_data_gau, num_bins = num_bins_vec[i], shape_val = shape_value)
  pred_df_training <- pred_df_training_object$pred_data
  centroid_df_training <- pred_df_training_object$df_bin_centroids
  avg_df_training <- pred_df_training_object$df_bin

  eval_df_training <- generate_eval_df(data = df_2, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, num_bins = num_bins_vec[i])

  pred_df_test_object <- predict_hex_id(training_data = training_data_5, nldr_df = tSNE_data_gau, nldr_df_test = tSNE_data_gau, num_bins = num_bins_vec[i], shape_val = shape_value)
  pred_df_test <- pred_df_test_object$pred_data
  centroid_df_test <- pred_df_test_object$df_bin_centroids
  avg_df_test <- pred_df_test_object$df_bin

  eval_df_test <- generate_eval_df(data = df_2, prediction_df = pred_df_test, df_bin_centroids = centroid_df_test, df_bin = avg_df_test, num_bins = num_bins_vec[i])

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)
  eval_data_test <- dplyr::bind_rows(eval_data_test, eval_df_test)

}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

eval_data_test <- eval_data_test |>
  mutate(data_type = "test")

MSE_df <- dplyr::bind_rows(eval_data_training, eval_data_test) |>
  dplyr::mutate(method = "tSNE")

## To draw with AIC
ggplot(MSE_df |> dplyr::filter(data_type == "training"), aes(x = number_of_bins,
                                                             y = total_error,
                                                             color = data_type
)) +
  geom_point() +
  geom_line() +
  #geom_vline(xintercept = NROW(full_grid_with_hexbin_id)) +
  #annotate("text", x= (NROW(full_grid_with_hexbin_id) - 10), y=-5000, label=paste0("effective number of bins = ", as.character(NROW(full_grid_with_hexbin_id))), angle=90) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02")) +
  ylab("AIC") +
  xlab("Total number of bins")
## Effective number of bins along x-axis

ggplot(MSE_df, aes(x = number_of_bins,
                   y = total_mse,
                   color = data_type
)) +
  geom_point() +
  geom_line() +
  # geom_vline(xintercept = NROW(full_grid_with_hexbin_id)) +
  # annotate("text", x= (NROW(full_grid_with_hexbin_id) - 10), y=0.25, label=paste0("effective number of bins = ", as.character(NROW(full_grid_with_hexbin_id))), angle=90) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02")) +
  ylab("MSE") +
  xlab("Total number of bins")



###########################################################

### UMAP

UMAP_fit <- umap(training_data_5 |> dplyr::select(-ID), n_neighbors = 15, n_components =  2)

UMAP_data_gau <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data_gau)[1:(ncol(UMAP_data_gau))] <- paste0(rep("UMAP",(ncol(UMAP_data_gau))), 1:(ncol(UMAP_data_gau)))

UMAP_data_gau <- UMAP_data_gau |>
  mutate(ID = training_data_5$ID)



#(perplexity: ", calculate_effective_perplexity(data), ")
plot_gau_umap <- plot_UMAP_2D(UMAP_data_gau) + #ggtitle("(a)") +
  theme_linedraw() +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())



cell_area <- 1

num_bins_umap_gau <- calculate_effective_x_bins(.data = UMAP_data_gau, x = UMAP1,
                                                cell_area = 1)
shape_val_umap_gau <- calculate_effective_shape_value(.data = UMAP_data_gau,
                                                      x = UMAP1, y = UMAP2)

## To extract bin centroids
hexbin_data_object_umap_gau <- extract_hexbin_centroids(nldr_df = UMAP_data_gau, num_bins = num_bins_umap_gau, shape_val = shape_val_umap_gau, x = UMAP1, y = UMAP2)

df_bin_centroids_umap_gau <- hexbin_data_object_umap_gau$hexdf_data

# min_std_cell_threshold3 <- 0.01
#
# df_bin_centroids_pbmc <- df_bin_centroids_pbmc |>
#   dplyr::mutate(stand_cell_count = counts/max(counts)) |>
#   dplyr::filter(stand_cell_count > min_std_cell_threshold3)

UMAP_data_with_hb_id_gau <- UMAP_data_gau |>
  dplyr::mutate(hb_id = hexbin_data_object_umap_gau$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_umap_gau <- dplyr::bind_cols(training_data_5 |> dplyr::select(-ID), UMAP_data_with_hb_id_gau)

## Averaged on high-D
df_bin_umap_gau <- avg_highD_data(.data = df_all_umap_gau)

## Triangulate bin centroids
tr1_object_umap_gau <- triangulate_bin_centroids(df_bin_centroids_umap_gau, x, y)
tr_from_to_df_umap_gau <- generate_edge_info(triangular_object = tr1_object_umap_gau)

## Compute 2D distances
distance_umap_gau <- cal_2D_dist(.data = tr_from_to_df_umap_gau)

## To find the benchmark value
benchmark_umap_gau <- find_benchmark_value(.data = distance_umap_gau, distance_col = distance)
#benchmark_tsne <- 5

tour_umap_gau <- show_langevitour(df_all_umap_gau, df_bin_umap_gau, df_bin_centroids_umap_gau, benchmark_value = benchmark_umap_gau, distance = distance_umap_gau, distance_col = distance)

## Prediction

shape_value <- calculate_effective_shape_value(.data = UMAP_data_gau,
                                               x = UMAP1, y = UMAP2)

## predict umap embeddings

predict_UMAP_df <- predict(UMAP_fit, test_data_5 |> dplyr::select(-ID)) |>
  as.data.frame()

names(predict_UMAP_df)[1:(ncol(predict_UMAP_df))] <- paste0(rep("UMAP",(ncol(predict_UMAP_df))), 1:(ncol(predict_UMAP_df)))

predict_UMAP_df <- predict_UMAP_df |>
  mutate(ID = test_data_5$ID)

plot_UMAP_2D(UMAP_data_gau) +
  geom_point(data = predict_UMAP_df, aes(x = UMAP1, y = UMAP2), color = "red")

num_bins_vec <- 1:20 ## Number of bins along the x-axis

vec <- stats::setNames(rep("", 6), c("number_of_bins", "number_of_observations", "total_error", "totol_error_method_2", "totol_error_method_3", "total_mse"))  ## Define column names

eval_data_test <- dplyr::bind_rows(vec)[0, ]
eval_data_test <- eval_data_test |>
  dplyr::mutate_if(is.character, as.numeric)

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(num_bins_vec)) {

  pred_df_training_object <- predict_hex_id(training_data = training_data_5, nldr_df = UMAP_data_gau, nldr_df_test = UMAP_data_gau, num_bins = num_bins_vec[i], shape_val = shape_value)
  pred_df_training <- pred_df_training_object$pred_data
  centroid_df_training <- pred_df_training_object$df_bin_centroids
  avg_df_training <- pred_df_training_object$df_bin

  eval_df_training <- generate_eval_df(data = df_2, prediction_df = pred_df_training, df_bin_centroids = centroid_df_training, df_bin = avg_df_training, num_bins = num_bins_vec[i])

  pred_df_test_object <- predict_hex_id(training_data = training_data_5, nldr_df = UMAP_data_gau, nldr_df_test = predict_UMAP_df, num_bins = num_bins_vec[i], shape_val = shape_value)
  pred_df_test <- pred_df_test_object$pred_data
  centroid_df_test <- pred_df_test_object$df_bin_centroids
  avg_df_test <- pred_df_test_object$df_bin

  eval_df_test <- generate_eval_df(data = df_2, prediction_df = pred_df_test, df_bin_centroids = centroid_df_test, df_bin = avg_df_test, num_bins = num_bins_vec[i])

  eval_data_training <- dplyr::bind_rows(eval_data_training, eval_df_training)
  eval_data_test <- dplyr::bind_rows(eval_data_test, eval_df_test)

}


## Add new column with data types

eval_data_training <- eval_data_training |>
  mutate(data_type = "training")

eval_data_test <- eval_data_test |>
  mutate(data_type = "test")

MSE_df <- dplyr::bind_rows(eval_data_training, eval_data_test) |>
  dplyr::mutate(method = "UMAP")

## To draw with AIC
ggplot(MSE_df |> dplyr::filter(data_type == "training"), aes(x = number_of_bins,
                                                             y = total_error,
                                                             color = data_type
)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = NROW(full_grid_with_hexbin_id)) +
  #annotate("text", x= (NROW(full_grid_with_hexbin_id) - 10), y=-5000, label=paste0("effective number of bins = ", as.character(NROW(full_grid_with_hexbin_id))), angle=90) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02")) +
  ylab("AIC") +
  xlab("Total number of bins")
## Effective number of bins along x-axis

ggplot(MSE_df, aes(x = number_of_bins,
                   y = total_mse,
                   color = data_type
)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = NROW(full_grid_with_hexbin_id)) +
  # annotate("text", x= (NROW(full_grid_with_hexbin_id) - 10), y=0.25, label=paste0("effective number of bins = ", as.character(NROW(full_grid_with_hexbin_id))), angle=90) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02")) +
  ylab("MSE") +
  xlab("Total number of bins")
