library(readr)
library(umap)
library(dplyr)

true_model <- read_rds("data/s_curve/true_model.rds") |>
  mutate(ID = row_number())

# ### UMAP
#
# UMAP_fit <- umap(true_model |> select(-ID), n_neighbors = 15, n_components =  2)
#
# UMAP_data_gau <- UMAP_fit$layout |>
#   as.data.frame()
# names(UMAP_data_gau)[1:(ncol(UMAP_data_gau))] <- paste0(rep("UMAP",(ncol(UMAP_data_gau))), 1:(ncol(UMAP_data_gau)))
#
# UMAP_data_gau <- UMAP_data_gau |>
#   mutate(ID = row_number())
#
# write_rds(UMAP_data_gau, file = "data/s_curve/umap_data_true_model.rds")
#
# scaled_umap <- gen_scaled_data(data = UMAP_data_gau)
#
# s_curve_noise_umap_scaled <- scaled_umap$scaled_nldr
# lim1 <- scaled_umap$lim1
# lim2 <- scaled_umap$lim2
# r2 <- diff(lim2)/diff(lim1)

# model <- fit_highd_model(training_data = true_model,
#                          emb_df = s_curve_noise_umap_scaled,
#                          bin1 = 12, r2 = r2,
#                          col_start_highd = "x")
# df_bin_centroids <- model$df_bin_centroids
# df_bin <- model$df_bin

pred_df_training <- predict_emb(test_data = true_model,
                             df_bin_centroids = df_bin_centroids_scurve,
                             df_bin = df_bin_scurve, type_NLDR = "UMAP")
glimpse(pred_df_training)

s_curve_noise_umap_scaled |>
    ggplot(aes(x = UMAP1,
               y = UMAP2,
               label = ID))+
    geom_point(alpha=0.5) +
    geom_point(data = pred_df_training, aes(x = pred_UMAP_1, y = pred_UMAP_2),
               color = "red", alpha=0.5) +
    coord_equal() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text = element_text(size = 5),
          axis.title = element_text(size = 7))


