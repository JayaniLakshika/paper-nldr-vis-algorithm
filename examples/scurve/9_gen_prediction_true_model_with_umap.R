library(readr)
library(umap)
library(dplyr)

true_model_df <- read_rds("data/s_curve/true_model.rds")
true_model_df <- true_model_df |>
  mutate(ID = row_number())

training_data_scurve <- read_rds("data/s_curve/s_curve_training.rds")

## To initialize effective bins along x
effective_bin1_scurve <- 12

scurve_model <- fit_highd_model(
  training_data = training_data_scurve,
  emb_df = umap_scurve_scaled,
  bin1 = effective_bin1_scurve,
  r2 = r2,
  q = 0.07,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_scurve <- scurve_model$df_bin_centroids
df_bin_scurve <- scurve_model$df_bin

true_pred_df <- predict_emb(
  test_data = true_model_df,
  df_bin_centroids = df_bin_centroids_scurve,
  df_bin = df_bin_scurve,
  type_NLDR = "UMAP"
)

true_pred_df |>
  ggplot(aes(x = pred_UMAP_1,
             y = pred_UMAP_2)) +
  geom_point(
    alpha=0.5,
    color = "#33a02c"
  )
