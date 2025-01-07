### This script is to generate the model withn tSNE for five Gaussian clusters data
library(readr)
library(dplyr)

training_data_gau <- read_rds("data/five_gau_clusters/data_five_gau.rds")

data_gau <- training_data_gau |>
  select(-ID) |>
  mutate(type = "data")

tsne_data_gau <- read_rds("data/five_gau_clusters/tsne_data_five_gau_71.rds")
gau1_scaled_obj <- gen_scaled_data(
  data = tsne_data_gau)
tsne_gau_scaled <- gau1_scaled_obj$scaled_nldr

tsne_gau <- tsne_gau_scaled |>
  ggplot(aes(x = tSNE1,
             y = tSNE2)) +
  geom_point(alpha=0.3)

## Compute hexbin parameters
num_bins_x_gau1 <- 20 #13
lim1 <- gau1_scaled_obj$lim1
lim2 <- gau1_scaled_obj$lim2
r2_gau1 <- diff(lim2)/diff(lim1)

gau1_model <- fit_highd_model(
  training_data = training_data_gau,
  emb_df = tsne_gau_scaled,
  bin1 = num_bins_x_gau1,
  r2 = r2_gau1,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x",
  q = 0.1
)

df_bin_centroids_gau1 <- gau1_model$df_bin_centroids
df_bin_gau1 <- gau1_model$df_bin

## Compute error

error_df_gau_abs <- augment(
  df_bin_centroids = df_bin_centroids_gau1,
  df_bin = df_bin_gau1,
  training_data = training_data_gau,
  newdata = NULL,
  type_NLDR = "tSNE",
  col_start = "x")

error_df_gau_abs <- error_df_gau_abs |>
  mutate(sqrt_row_wise_abs_error = sqrt(row_wise_abs_error))

error_df_gau_abs <- error_df_gau_abs |>
  bind_cols(tsne_gau_scaled |>
              select(-ID))

error_df_gau_tsne <- error_df_gau_abs |>
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             z = sqrt_row_wise_abs_error)) +
  stat_summary_hex(bins=15) +
  xlim(c(0.3, 0.66)) + ylim(c(-0.01, 0.35)) +
  scale_fill_continuous_sequential(palette = "YlOrRd", n_interp = 20) +
  interior_annotation("a3",
                      position = c(0.08, 0.9),
                      cex = 2)

## Triangulate bin centroids
tr1_object_gau1 <- tri_bin_centroids(
  df_bin_centroids_gau1, x = "c_x", y = "c_y")
tr_from_to_df_gau1 <- gen_edges(
  tri_object = tr1_object_gau1)

# tr_from_to_df_gau1 <- tr_from_to_df_gau1 |>
#   filter(row_number() != 76) |>
#   filter(row_number() != 34)

# tr_from_to_df_gau1 <- tr_from_to_df_gau1 |>
#   filter(row_number() != 149)

## Compute 2D distances
distance_gau1 <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_gau1,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_gau1 <- find_lg_benchmark(
  distance_edges = distance_gau1,
  distance_col = "distance")

benchmark_gau1 <- 0.1

# trimesh_removed_gau1 <- vis_rmlg_mesh(
#   distance_edges = distance_gau1,
#   benchmark_value = benchmark_gau1,
#   tr_coord_df = tr_from_to_df_gau1,
#   distance_col = "distance")

tr_df <- distinct(tibble(
  x = c(tr_from_to_df_gau1[["x_from"]], tr_from_to_df_gau1[["x_to"]]),
  y = c(tr_from_to_df_gau1[["y_from"]], tr_from_to_df_gau1[["y_to"]])))

distance_df_small_edges_gau1 <- distance_gau1 |>
  filter(distance < benchmark_gau1) |>
  filter(!(row_number() %in% c(33, 40, 48, 138, 143)))

tr_from_to_df_gau1 <- inner_join(
  tr_from_to_df_gau1, distance_df_small_edges_gau1,
  by = c("from", "to"))

trimesh_removed_gau_tsne <- ggplot() +
  geom_point(
    data = tsne_gau_scaled,
    aes(
      x = tSNE1,
      y = tSNE2
    ),
    alpha = 0.2) +
  geom_rect(
    aes(xmin = 0.3,
        xmax = 0.65,
        ymin = -0.05,
        ymax = 0.35),
    fill = "transparent",
    color = "grey50",
    linewidth = 0.5) +
  interior_annotation("a1",
                      position = c(0.08, 0.95),
                      cex = 2) +
  theme(
    aspect.ratio = 1
  )


## Hexagonal binning to have regular hexagons
hb_obj_gau1 <- hex_binning(
  data = tsne_gau_scaled,
  bin1 = num_bins_x_gau1,
  r2 = r2_gau1,
  q = 0.1)

tsne_data_with_hb_id <- hb_obj_gau1$data_hb_id

df_all_gau1 <- dplyr::bind_cols(training_data_gau |> dplyr::select(-ID),
                                tsne_data_with_hb_id)

### Define type column
df <- df_all_gau1 |>
  dplyr::select(tidyselect::starts_with("x")) |>
  dplyr::mutate(type = "data") ## original dataset

df_b <- df_bin_gau1 |>
  dplyr::filter(hb_id %in% df_bin_centroids_gau1$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_gau1$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

# Apply the scaling
df_model_data <- bind_rows(data_gau, df_b)
scaled_gau <- scale_data_manual(df_model_data, "type") |>
  as_tibble()

scaled_gau_data <- scaled_gau |>
  filter(type == "data") |>
  select(-type)

scaled_gau_data_model <- scaled_gau |>
  filter(type == "model") |>
  select(-type)

## Set the maximum difference as the criteria
distance_df_small_edges <- distance_gau1 |>
  dplyr::filter(distance < benchmark_gau1)

projection <- cbind(
  c(-0.00215,-0.68905,-0.04778,-0.54223),
  c(0.42558,-0.23854,-0.63659,0.35753))

projection_scaled <- projection * 1

projected <- as.matrix(scaled_gau_data) %*% projection_scaled

projected_df <- projected |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::rename(c("proj1" = "...1",
                  "proj2" = "...2")) |>
  #dplyr::mutate(type = df_exe$type) |>
  dplyr::mutate(ID = dplyr::row_number())

projected_model <- as.matrix(scaled_gau_data_model) %*% projection_scaled

projected_model_df <- projected_model |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::rename(c("proj1" = "...1",
                  "proj2" = "...2")) |>
  dplyr::mutate(ID = dplyr::row_number())

model_df <- dplyr::left_join(
  distance_df_small_edges |> select(-distance),
  projected_model_df,
  by = c("from" = "ID"))

names(model_df)[3:NCOL(model_df)] <- paste0(names(projected_model_df)[-NCOL(projected_model_df)], "_from")

model_df <- dplyr::left_join(model_df, projected_model_df, by = c("to" = "ID"))
names(model_df)[(2 + NCOL(projected_model_df)):NCOL(model_df)] <- paste0(names(projected_model_df)[-NCOL(projected_model_df)], "_to")

axes_obj <- gen_axes(
  proj = projection,
  limits = 0.3,
  axis_pos_x = -0.6,
  axis_pos_y = -0.6,
  axis_labels = names(scaled_gau_data),
  threshold = 0)

axes <- axes_obj$axes
circle <- axes_obj$circle

five_gau_proj_tsne_model1 <- projected_df |>
  ggplot(
    aes(
      x = proj1,
      y = proj2)) +
  geom_segment(
    data = model_df,
    aes(
      x = proj1_from,
      y = proj2_from,
      xend = proj1_to,
      yend = proj2_to),
    color = "#000000",
    linewidth = 1) +
  geom_point(
    #size = 0.5,
    alpha = 0.2,
    color = "#999999") +
  geom_segment(
    data=axes,
    aes(x=x1, y=y1, xend=x2, yend=y2),
    colour="grey70") +
  geom_text(
    data=axes,
    aes(x=x2, y=y2),
    label=rownames(axes),
    colour="grey50",
    size = 5) +
  geom_path(
    data=circle,
    aes(x=c1, y=c2), colour="grey70") +
  coord_fixed() +
  xlim(c(-0.83, 0.83)) +
  ylim(c(-0.83, 0.83)) +
  interior_annotation("a2",
                      position = c(0.08, 0.9),
                      cex = 2)

gen_proj_five_clust <- function(){
  list(tsne_gau, five_gau_proj_tsne_model1)
}
