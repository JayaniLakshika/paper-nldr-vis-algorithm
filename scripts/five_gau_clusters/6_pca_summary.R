## This script is to check PCA results for Five Gaussian clusters data
library(tidyverse)
library(quollr)

calculate_pca <- function(feature_dataset){
  pcaY_cal <- prcomp(feature_dataset, center = TRUE, scale = FALSE)
  PCAresults <- data.frame(pcaY_cal$x[, 1:4])
  summary_pca <- summary(pcaY_cal)
  var_explained_df <- data.frame(PC= paste0("PC",1:4),
                                 var_explained=(pcaY_cal$sdev[1:2])^2/sum((pcaY_cal$sdev[1:2])^2))
  return(list(prcomp_out = pcaY_cal,pca_components = PCAresults, summary = summary_pca, var_explained_pca  = var_explained_df))
}

data <- read_rds("data/five_gau_clusters/data_five_gau_with_clusts.rds")
data_n <- data
data <- data |>
  dplyr::select(-ID)

pca_ref_calc <- calculate_pca(data |> dplyr::select(where(is.numeric)))
data_pca <- pca_ref_calc$pca_components |>
  dplyr::mutate(ID = row_number()) |>
  dplyr::mutate(cluster = data$cluster)

# gau1_scaled_obj <- gen_scaled_data(
#   data = data_pca |> dplyr::select(-cluster))
# pca_gau_scaled <- gau1_scaled_obj$scaled_nldr
#
# ## Compute hexbin parameters
# num_bins_x_gau1 <- 14
# lim1 <- gau1_scaled_obj$lim1
# lim2 <- gau1_scaled_obj$lim2
# r2_gau1 <- diff(lim2)/diff(lim1)
#
# gau1_model <- fit_highd_model(
#   highd_data = data_n,
#   nldr_data = pca_gau_scaled,
#   bin1 = num_bins_x_gau1,
#   r2 = r2_gau1,
#   is_bin_centroid = TRUE,
#   q = 0.1
# )
#
# df_bin_centroids_gau1 <- gau1_model$df_bin_centroids
# df_bin_gau1 <- gau1_model$df_bin
#
# ## Triangulate bin centroids
# tr1_object_gau1 <- tri_bin_centroids(
#   df_bin_centroids_gau1, x = "c_x", y = "c_y")
# tr_from_to_df_gau1 <- gen_edges(
#   tri_object = tr1_object_gau1)
#
# ## Compute 2D distances
# distance_gau1 <- cal_2d_dist(
#   tr_coord_df = tr_from_to_df_gau1,
#   start_x = "x_from",
#   start_y = "y_from",
#   end_x = "x_to",
#   end_y = "y_to",
#   select_vars = c("from", "to", "distance"))
#
# ## To find the benchmark value
# benchmark_gau1 <- find_lg_benchmark(
#   distance_edges = distance_gau1,
#   distance_col = "distance")
#
# tr_df <- distinct(tibble(
#   x = c(tr_from_to_df_gau1[["x_from"]], tr_from_to_df_gau1[["x_to"]]),
#   y = c(tr_from_to_df_gau1[["y_from"]], tr_from_to_df_gau1[["y_to"]])))
#
# distance_df_small_edges_gau1 <- distance_gau1 |>
#   filter(distance < 0.1) #benchmark_gau1
#
# distance_df_small_edges_gau1 <- distance_df_small_edges_gau1 |>
#   mutate(ID = row_number()) |>
#   dplyr::filter(!(ID %in% c(150, 110))) |>
#   dplyr::select(-ID)
#
# tr_from_to_df_gau1 <- right_join(
#   tr_from_to_df_gau1, distance_df_small_edges_gau1,
#   by = c("from", "to"))
#
# ## Hexagonal binning to have regular hexagons
# hb_obj_gau1 <- hex_binning(
#   data = tsne_gau_scaled,
#   bin1 = num_bins_x_gau1,
#   r2 = r2_gau1,
#   q = 0.1)
#
# tsne_data_with_hb_id <- hb_obj_gau1$data_hb_id
#
# df_all_gau1 <- dplyr::bind_cols(training_data_gau |> dplyr::select(-ID),
#                                 tsne_data_with_hb_id)
#
# ### Define type column
# df <- df_all_gau1 |>
#   dplyr::select(tidyselect::starts_with("x")) |>
#   dplyr::mutate(type = "data") ## original dataset
#
# df_b <- df_bin_gau1 |>
#   dplyr::filter(hb_id %in% df_bin_centroids_gau1$hexID) |>
#   dplyr::mutate(type = "model") ## Data with summarized mean
#
# ## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
# df_b <- df_b[match(df_bin_centroids_gau1$hexID, df_b$hb_id),] |>
#   dplyr::select(-hb_id)
#
# # Apply the scaling
# df_model_data <- bind_rows(data_gau, df_b)
# scaled_gau <- scale_data_manual(df_model_data, "type") |>
#   as_tibble()
#
# scaled_gau_data <- scaled_gau |>
#   filter(type == "data") |>
#   select(-type)
#
# scaled_gau_data_model <- scaled_gau |>
#   filter(type == "model") |>
#   select(-type)



## PC1 Vs PC2

ggplot(data_pca, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5) +
  theme(
    aspect.ratio = 1
  )

## PC1 Vs PC3

ggplot(data_pca, aes(x = PC1, y = PC3)) +
  geom_point(alpha = 0.5) +
  theme(
    aspect.ratio = 1
  )

## PC3 Vs PC4

ggplot(data_pca, aes(x = PC3, y = PC4)) +
  geom_point(alpha = 0.5) +
  theme(
    aspect.ratio = 1
  )

## For selected cluster

data_pca_cluster1 <- data_pca |>
  dplyr::filter(cluster == "cluster1")

ggplot(data_pca_cluster1, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5) +
  theme(
    aspect.ratio = 1
  )

## PC1 Vs PC3

ggplot(data_pca_cluster1, aes(x = PC1, y = PC3)) +
  geom_point(alpha = 0.5) +
  theme(
    aspect.ratio = 1
  )

## PC3 Vs PC4

ggplot(data_pca_cluster1, aes(x = PC3, y = PC4)) +
  geom_point(alpha = 0.5) +
  theme(
    aspect.ratio = 1
  )
