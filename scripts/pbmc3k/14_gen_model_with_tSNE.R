## This script is to generate the model with tSNE for PBMC3k data

library(quollr)
library(tidyverse)

source("scripts/additional_functions.R")
set.seed(20240110)

## Data

training_data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50.rds")
names(training_data_pbmc) <- paste0("x", 1:50)

training_data_pbmc <- training_data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(training_data_pbmc))

## NLDR
tsne_pbmc <- read_rds("data/pbmc3k/pbmc_tsne_30.rds")

## To fit the model

## Compute hexbin parameters
num_bins_x_tsne_pbmc <- 22

algo_obj_pbmc <- fit_highd_model(
  highd_data = training_data_pbmc,
  nldr_data = tsne_pbmc,
  b1 = num_bins_x_tsne_pbmc,
  q = 0.1,
  benchmark_highdens = 1)

tsne_pbmc_scaled_best <- algo_obj_pbmc$nldr_obj$scaled_nldr
tr_from_to_df_pbmc <- algo_obj_pbmc$trimesh_data
df_bin_centroids_pbmc <- algo_obj_pbmc$model_2d
df_bin_pbmc <- algo_obj_pbmc$model_highd

tsne_pbmc_scaled_best_with_cluster <- tsne_pbmc_scaled_best |>
  mutate(cluster = umap_pbmc_scaled_with_cluster$cluster)

write_rds(tsne_pbmc_scaled_best_with_cluster, "data/pbmc3k/tsne_pbmc_scaled_best_with_cluster.rds")
write_rds(tr_from_to_df_pbmc, "data/pbmc3k/tsne_tr_from_to_df_pbmc.rds")

df_b_pbmc <- df_bin_pbmc |>
  dplyr::filter(h %in% df_bin_centroids_pbmc$h) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the h order in df_b_with_center_data
df_b_pbmc <- df_b_pbmc[match(df_bin_centroids_pbmc$h, df_b_pbmc$h),] |>
  dplyr::select(-h)

# Apply the scaling
df_model_data_pbmc <- bind_rows(data_pbmc, df_b_pbmc)
scaled_pbmc <- scale_data_manual(df_model_data_pbmc, "type") |>
  as_tibble()

scaled_pbmc_data <- scaled_pbmc |>
  filter(type == "data") |>
  select(-type)

scaled_pbmc_data_model <- scaled_pbmc |>
  filter(type == "model") |>
  select(-type)

df_model_data_pbmc_n <- bind_rows(df_b_pbmc, data_pbmc_n)

langevitour::langevitour(df_model_data_pbmc_n[1:(length(df_model_data_pbmc_n)-1)],
                         lineFrom = tr_from_to_df_pbmc$from,
                         lineTo = tr_from_to_df_pbmc$to,
                         group = factor(df_model_data_pbmc_n$type,
                                        c("cluster1", "cluster2", "cluster3", "model")),
                         levelColors = c("#8dd3c7", "#fdb462", "#bebada", "#000000"))

## Model projections

## First projection
projection <- cbind(
  c(-0.4545,-0.5533,-0.1235,0.3984,-0.2065,-0.3586,-0.0088,0.3572,-0.1368),
  c(0.4576,-0.3157,0.4708,0.1410,-0.4479,-0.1685,0.3668,-0.2692,0.1332))

proj_obj1 <- get_projection(projection = projection,
                            proj_scale = 1.23,
                            highd_data = scaled_pbmc_data,
                            model_highd = scaled_pbmc_data_model,
                            trimesh_data = tr_from_to_df_pbmc,
                            axis_param = list(limits = 1,
                                              axis_scaled = 1,
                                              axis_pos_x = -0.48,
                                              axis_pos_y = -0.48,
                                              threshold = 0.082))

proj_obj1[["cluster"]] <- as.character(umap_pbmc_scaled_with_cluster$cluster)

write_rds(proj_obj1, "data/pbmc3k/pbmc_tsne_proj_obj1.rds")

## Second projection
projection <- cbind(
  c(0.0734,0.6138,-0.4671,0.2589,0.1189,0.0778,-0.0819,0.0308,-0.5561),
  c(0.6625,0.0034,0.2412,0.3360,0.2463,0.1070,-0.2903,-0.4705,0.1293))

proj_obj2 <- get_projection(projection = projection,
                            proj_scale = 1.25,
                            highd_data = scaled_pbmc_data,
                            model_highd = scaled_pbmc_data_model,
                            trimesh_data = tr_from_to_df_pbmc,
                            axis_param = list(limits = 1.35,
                                              axis_scaled = 1,
                                              axis_pos_x = -0.75,
                                              axis_pos_y = -0.75,
                                              threshold = 0.12))

proj_obj2[["cluster"]] <- as.character(umap_pbmc_scaled_with_cluster$cluster)

write_rds(proj_obj2, "data/pbmc3k/pbmc_tsne_proj_obj2.rds")

