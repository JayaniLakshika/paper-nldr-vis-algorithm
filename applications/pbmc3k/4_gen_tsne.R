library(readr)
library(Rtsne)
library(dplyr)
library(ggplot2)
library(quollr)
set.seed(20240110)

perplexity <- 51

## Select PCs
training_data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50.rds")
training_data_pbmc <- training_data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(training_data_pbmc))

## Obtain tSNE
tSNE_data <- Fit_tSNE(training_data_pbmc |> select(-ID), opt_perplexity = perplexity, with_seed = 20240110)

tSNE_data <- tSNE_data |>
  select(-ID) |>
  mutate(ID = training_data_pbmc$ID)

## Run only once
write_rds(tSNE_data, file = paste0("data/pbmc3k/pbmc_tsne_", perplexity, ".rds"))

## Scaled data
tsne_pbmc_scaled <- as.data.frame(do.call(cbind, gen_scaled_data(data = tSNE_data,
                                                                 x = "tSNE1", y = "tSNE2"))) |>
  dplyr::rename(c("tSNE1" = "scaled_tSNE1",
                  "tSNE2" = "scaled_tSNE2")) |>
  dplyr::mutate(ID = 1:NROW(tSNE_data))

## UMAP
## Prediction

hex_size_vec <- seq(0.01, 2, by = 0.01)

vec <- stats::setNames(rep("", 7), c("num_bins", "error", "mse", "num_bins_x", "num_bins_y", "hex_size", "num_non_empty_bins"))  ## Define column names

mse_df_pbmc <- dplyr::bind_rows(vec)[0, ]
mse_df_pbmc <- mse_df_pbmc |>
  dplyr::mutate_if(is.character, as.numeric)

for (i in 1:length(hex_size_vec)) {

  num_bin_list <- calc_bins(data = tsne_pbmc_scaled,
                            x = "tSNE1", y = "tSNE2",
                            hex_size = hex_size_vec[i], buffer_x = NA, buffer_y = NA)

  num_bins_x <- num_bin_list$num_x
  num_bins_y <- num_bin_list$num_y

  model_object <- fit_highd_model( training_data = training_data_pbmc,
                                   nldr_df_with_id = tsne_pbmc_scaled,
                                   x = "tSNE1", y = "tSNE2",
                                   num_bins_x = num_bins_x,
                                   num_bins_y = num_bins_y,
                                   x_start = NA, y_start = NA,
                                   buffer_x = NA, buffer_y = NA,
                                   hex_size = hex_size_vec[i],
                                   is_rm_lwd_hex = FALSE,
                                   benchmark_to_rm_lwd_hex = NA,
                                   col_start_2d = "tSNE",
                                   col_start_highd = "PC_")

  centroid_df_training <- model_object$df_bin_centroids
  avg_df_training <- model_object$df_bin

  pred_emb_list <- predict_emb(test_data = training_data_pbmc,
                               df_bin_centroids = centroid_df_training,
                               df_bin = avg_df_training, type_NLDR = "tSNE")

  pred_df_training <- as.data.frame(do.call(cbind, pred_emb_list))

  eval_list <- gen_summary(test_data = training_data_pbmc,
                           prediction_df = pred_df_training,
                           df_bin = avg_df_training, col_start = "PC_")

  mse_df_pbmc <- mse_df_pbmc |>
    tibble::add_row(num_bins = num_bins_x * num_bins_y,
                    error = eval_list$error,
                    mse = eval_list$mse,
                    num_bins_x = num_bins_x,
                    num_bins_y = num_bins_y,
                    hex_size = hex_size_vec[i],
                    num_non_empty_bins = NROW(centroid_df_training))


}


## If same total number of bins occurred only select ones with minimum error
### Obtain duplicate bins
dupli_bins <- mse_df_pbmc |>
  dplyr::count(num_bins) |>
  dplyr::filter(n > 1) |>
  dplyr::pull(num_bins)

### Group split by duplicated bins
duplicate_df_list <- mse_df_pbmc |>
  dplyr::filter(num_bins %in% dupli_bins) |>
  dplyr::arrange(num_bins) |>
  dplyr::group_split(num_bins)

### Obtain one row from duplicates which have lowest error and hexsize
duplicate_df <- data.frame(matrix(nrow = 0, ncol = 0))

for (i in 1:length(duplicate_df_list)) {

  dd <- duplicate_df_list[[i]] |>
    dplyr::filter(mse == min(duplicate_df_list[[i]]$mse)) |>
    dplyr::filter(hex_size == min(duplicate_df_list[[i]]$hex_size))

  duplicate_df <- dplyr::bind_rows(duplicate_df, dd)

}

### Obtain the mse_df with not duplicated bins
not_dupli_df <- mse_df_pbmc |>
  dplyr::filter(!(num_bins %in% dupli_bins))

### Combine duplicated and not duplicated(corrected) bins dfs
mse_df_pbmc <- dplyr::bind_rows(not_dupli_df, duplicate_df) |>
  dplyr::mutate(perplexity = perplexity)


abs_error_plot_pbmc <- ggplot(mse_df_pbmc, aes(x = num_non_empty_bins,
                                               y = error
)) +
  geom_point() +
  geom_line() +
  # geom_vline(xintercept = 588, linetype="solid",
  #              color = "red", size=0.8, alpha = 0.5) +
  theme_light() +
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7)) +
  ylab("absolute error") +
  xlab("number of non-empty bins")

abs_error_plot_pbmc

write_csv(mse_df_pbmc, "data/pbmc3k/error_df_tsne.csv", append = TRUE)

