library(tidyverse)
library(quollr)

tsne_pbmc <- read_rds("data/pbmc3k/pbmc_scdeed_tsne_perplexity_30.rds")
tsne_pbmc <- as_tibble(tsne_pbmc)
names(tsne_pbmc) <- c("emb1", "emb2")
tsne_pbmc <- tsne_pbmc |>
  dplyr::mutate(ID = row_number())


## Compute hexbin parameters
num_bins_x_pbmc <- 30

algo_obj_pbmc <- fit_highd_model(
  highd_data = training_data_pbmc,
  nldr_data = tsne_pbmc,
  b1 = num_bins_x_pbmc,
  q = 0.1,
  benchmark_highdens = 5)

tsne_pbmc_scaled <- algo_obj_pbmc$nldr_obj$scaled_nldr
#tr_from_to_df_pbmc <- algo_obj_pbmc$trimesh_data
df_bin_centroids_pbmc <- algo_obj_pbmc$model_2d
df_bin_pbmc <- algo_obj_pbmc$model_highd

## Removed some bins manually to fit a better model
df_bin_centroids_pbmc <- df_bin_centroids_pbmc |>
  filter(!(hexID %in% c(559, 617, 677)))

tri_object <- tri_bin_centroids(
  centroids_data = df_bin_centroids_pbmc)

tr_from_to_df_pbmc <- gen_edges(tri_object = tri_object, a1 = algo_obj_pbmc$hb_obj$a1)
tr_from_to_df_pbmc <- update_trimesh_index(trimesh_data = tr_from_to_df_pbmc)


tsne_pbmc_scaled_with_cluster <- tsne_pbmc_scaled |>
  mutate(cluster = if_else((emb2 <= -emb1 + 0.8) & (emb1 < 0.5), "cluster1", if_else(emb2 >= -2 * emb1 + 1.79, "cluster3", "cluster2"))) |>
  mutate(cluster = factor(cluster, levels = c("cluster1", "cluster2", "cluster3")))

write_rds(tsne_pbmc_scaled_with_cluster, "data/pbmc3k/tsne_pbmc_scaled_with_cluster_scDEED.rds")
write_rds(tr_from_to_df_pbmc, "data/pbmc3k/tr_from_to_df_pbmc_tsne_scDEED.rds")

## 2D projection

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
  mutate(type = as.character(tsne_pbmc_scaled_with_cluster$cluster))

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

proj_obj1[["cluster"]] <- as.character(tsne_pbmc_scaled_with_cluster$cluster)

write_rds(proj_obj1, "data/pbmc3k/proj_obj1_tsne_scDEED.rds")


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

proj_obj2[["cluster"]] <- as.character(tsne_pbmc_scaled_with_cluster$cluster)

write_rds(proj_obj2, "data/pbmc3k/proj_obj2_tsne_scDEED.rds")
