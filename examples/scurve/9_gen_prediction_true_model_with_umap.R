# Load necessary libraries
library(readr)
library(umap)
library(dplyr)
library(ggplot2)

# Read data
true_model_df <- read_rds("data/s_curve/scurve_true_model.rds")
training_data_scurve <- read_rds("data/s_curve/s_curve_training.rds")

# Initialize effective bins along x
effective_bin1_scurve <- 15

# Fit high-dimensional model
scurve_model <- fit_highd_model(
  training_data = training_data_scurve,
  emb_df = umap_scurve_scaled,
  bin1 = effective_bin1_scurve,
  r2 = r2,
  q = 0.16,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_scurve <- scurve_model$df_bin_centroids
df_bin_scurve <- scurve_model$df_bin

# Predict embeddings
true_pred_df <- predict_emb(
  test_data = true_model_df,
  df_bin_centroids = df_bin_centroids_scurve,
  df_bin = df_bin_scurve,
  type_NLDR = "UMAP"
)

# Compute hexbin parameters
num_bins_x_scurve <- 15

# Hexagon binning to have regular hexagons
hb_obj_scurve <- hex_binning(
  data = umap_scurve_scaled,
  bin1 = num_bins_x_scurve,
  r2 = r2,
  q = 0.16
)

# Data set with all centroids
all_centroids_df <- hb_obj_scurve$centroids

# Generate all coordinates of hexagons
hex_grid <- hb_obj_scurve$hex_poly

# Obtain the standardized counts within hexbins
counts_df <- hb_obj_scurve$std_cts

hex_grid_with_counts <- left_join(hex_grid, counts_df, by = c("hex_poly_id" = "hb_id"))
hex_grid_nonempty <- hex_grid %>% filter(hex_poly_id %in% df_bin_centroids_scurve$hexID)

# Compute radius r
a1 <- calc_bins_y(bin1 = 15, r2 = r2, q = 0.16)$a1
r <- a1 / 2

# Add jitter to the predicted points within each hexagon
set.seed(123)  # for reproducibility

# Function to jitter points within a circumcircle of a hexagon
jitter_within_circumcircle <- function(center_x, center_y, radius, num_points) {
  theta <- runif(num_points, 0, 2 * pi)
  r <- radius * sqrt(runif(num_points))
  x <- center_x + r * cos(theta)
  y <- center_y + r * sin(theta)
  jittered_points <- data.frame(x = x, y = y)
  return(jittered_points)
}

# Jittered points data frame
jittered_points_df <- data.frame()

# Iterate through each hexagon and jitter points within the circumcircle
for (hex_id in unique(hex_grid_nonempty$hex_poly_id)) {
  hex_points <- hex_grid_nonempty %>% filter(hex_poly_id == hex_id)
  center_x <- mean(hex_points$x)
  center_y <- mean(hex_points$y)
  points_in_hex <- true_pred_df %>% filter(pred_UMAP_1 >= min(hex_points$x) & pred_UMAP_1 <= max(hex_points$x) & pred_UMAP_2 >= min(hex_points$y) & pred_UMAP_2 <= max(hex_points$y))

  if (nrow(points_in_hex) > 0) {
    jittered_points <- jitter_within_circumcircle(center_x, center_y, r, nrow(points_in_hex))
    jittered_points_df <- rbind(jittered_points_df, jittered_points)
  }
}

# Plotting with jittered points within hexagons
ggplot(data = hex_grid_nonempty, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = hex_poly_id), fill = "#ffffff") +
  geom_point(data = jittered_points_df, aes(x = x, y = y), alpha = 0.5, color = "#969696")
