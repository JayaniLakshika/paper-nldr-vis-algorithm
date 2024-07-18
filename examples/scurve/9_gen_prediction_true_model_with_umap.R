library(readr)
library(umap)
library(dplyr)

true_model_df <- read_rds("data/s_curve/true_model.rds")
true_model_df <- true_model_df |>
  mutate(ID = row_number())

training_data_scurve <- read_rds("data/s_curve/s_curve_training.rds")

## To initialize effective bins along x
effective_bin1_scurve <- 11

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

## Compute hexbin parameters
num_bins_x_scurve <- 11

## hexagon binning to have regular hexagons
hb_obj_scurve <- hex_binning(
  data = umap_scurve_scaled,
  bin1 = num_bins_x_scurve,
  r2 = r2,
  q = 0.07)

## Data set with all centroids
all_centroids_df <- hb_obj_scurve$centroids

## Generate all coordinates of hexagons
hex_grid <- hb_obj_scurve$hex_poly

## To obtain the standardise counts within hexbins
counts_df <- hb_obj_scurve$std_cts

hex_grid_with_counts <-
  left_join(hex_grid,
            counts_df,
            by = c("hex_poly_id" = "hb_id"))

hex_grid_nonempty <- hex_grid |>
  filter(hex_poly_id %in% df_bin_centroids_scurve$hexID)

ggplot(
  data = hex_grid_nonempty,
  aes(x = x, y = y)) +
  geom_polygon(color = "black",
               aes(group = hex_poly_id),
               fill = "#ffffff") +
  geom_jitter(
    data = true_pred_df,
    aes(x = pred_UMAP_1,
        y = pred_UMAP_2),
    alpha=0.5,
    color = "#969696"
  )
