library(tidyverse)
library(quollr)

## NLDR data
# result_umap <- read_rds(here::here('data/pbmc3k/pbmc_scdeed_umap_results.rds'))
#
# result_umap |>
#   filter(min.dist == 0.5) |>
#   filter(n_neighbors < 240) |>
#   ggplot(aes(
#     x = as.numeric(n_neighbors),
#     y = number_dubious_cells
#   )) +
#   geom_line() +
#   geom_point() +
#   labs(x = "n_neighbors",
#        y = "Number of dubious cells") +
#   theme_minimal()

umap_pbmc <- read_rds("data/pbmc3k/pbmc_scdeed_umap_n_neighbors_30_min_dist_0.3.rds")
umap_pbmc <- as_tibble(umap_pbmc)
names(umap_pbmc) <- c("emb1", "emb2")
umap_pbmc <- umap_pbmc |>
  dplyr::mutate(ID = row_number())

training_data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50_scdeed.rds")
training_data_pbmc <- training_data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(training_data_pbmc))

names(training_data_pbmc) <- append(paste0("x", 1:9), "ID")

## Compute hexbin parameters
num_bins_x_pbmc <- 30

algo_obj_pbmc <- fit_highd_model(
  highd_data = training_data_pbmc,
  nldr_data = umap_pbmc,
  bin1 = num_bins_x_pbmc,
  q = 0.1,
  benchmark_highdens = 5)

umap_pbmc_scaled <- algo_obj_pbmc$nldr_obj$scaled_nldr
tr_from_to_df_pbmc <- algo_obj_pbmc$trimesh_data
df_bin_centroids_pbmc <- algo_obj_pbmc$model_2d
df_bin_pbmc <- algo_obj_pbmc$model_highd

umap_pbmc_scaled_with_cluster <- umap_pbmc_scaled |>
  mutate(cluster = if_else(emb1 <= 0.3, "cluster1", if_else((emb1 >= 0.75) & (emb2 <= 0.25), "cluster3", "cluster2"))) |>
  mutate(cluster = factor(cluster, levels = c("cluster1", "cluster2", "cluster3")))

write_rds(umap_pbmc_scaled_with_cluster, "data/pbmc3k/umap_pbmc_scaled_with_cluster_scDEED.rds")
write_rds(tr_from_to_df_pbmc, "data/pbmc3k/tr_from_to_df_pbmc_scDEED.rds")


## 2D projections

data_pbmc <- training_data_pbmc |>
  select(-ID) |>
  mutate(type = "data")

df_b_pbmc <- df_bin_pbmc |>
  dplyr::filter(hexID %in% df_bin_centroids_pbmc$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b_pbmc <- df_b_pbmc[match(df_bin_centroids_pbmc$hexID, df_b_pbmc$hexID),] |>
  dplyr::select(-hexID)

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

data_pbmc_n <- data_pbmc |>
  select(-type) |>
  mutate(type = as.character(umap_pbmc_scaled_with_cluster$cluster))

df_model_data_pbmc_n <- bind_rows(df_b_pbmc, data_pbmc_n)

langevitour::langevitour(df_model_data_pbmc_n[1:(length(df_model_data_pbmc_n)-1)],
                         lineFrom = tr_from_to_df_pbmc$from,
                         lineTo = tr_from_to_df_pbmc$to,
                         group = factor(df_model_data_pbmc_n$type,
                                        c("cluster1", "cluster2", "cluster3", "model")),
                         levelColors = c("#8dd3c7", "#fdb462", "#bebada", "#000000"))

## First projection
projection <- cbind(
  c(0.025588,0.002516,0.002185,0.025284,0.004030,-0.026880,-0.015326,0.014851,0.026260),
  c(0.026495,-0.024817,0.021623,-0.014036,-0.010220,0.016634,-0.006496,-0.022928,0.016045))

proj_obj1 <- get_projection(projection = projection,
                            proj_scale = 1.23,
                            highd_data = scaled_pbmc_data,
                            model_highd = scaled_pbmc_data_model,
                            trimesh_data = tr_from_to_df_pbmc,
                            axis_param = list(limits = 0.05,
                                              axis_scaled = 17,
                                              axis_pos_x = -0.03,
                                              axis_pos_y = -0.03,
                                              threshold = 0.0044))

proj_obj1[["cluster"]] <- as.character(umap_pbmc_scaled_with_cluster$cluster)

write_rds(proj_obj1, "data/pbmc3k/proj_obj1_umap_scDEED.rds")

## Second projection
projection <- cbind(
  c(0.025013,0.029055,-0.003143,0.001821,-0.006208,0.019174,0.020630,0.021512,0.020281),
  c(0.004528,0.006465,-0.042452,0.009837,0.010040,-0.023820,0.005236,0.014775,-0.017712))

proj_obj2 <- get_projection(projection = projection,
                            proj_scale = 1.25,
                            highd_data = scaled_pbmc_data,
                            model_highd = scaled_pbmc_data_model,
                            trimesh_data = tr_from_to_df_pbmc,
                            axis_param = list(limits = 0.035,
                                              axis_scaled = 17,
                                              axis_pos_x = -0.023,
                                              axis_pos_y = -0.023,
                                              threshold = 0.0027))

proj_obj2[["cluster"]] <- as.character(umap_pbmc_scaled_with_cluster$cluster)

write_rds(proj_obj2, "data/pbmc3k/proj_obj2_umap_scDEED.rds")
