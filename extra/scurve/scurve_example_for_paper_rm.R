```{r}
#| label: scurve-training
training_data_scurve <- read_rds("data/s_curve/s_curve_training.rds")
```
```{r}
#| label: scurve-training
training_data_scurve <- read_rds("data/s_curve/s_curve_training.rds")
```

<!--UMAP applied for Scurve data-->
  ```{r}
#| label: umap-scurve
umap_scurve <- read_rds(file = "data/s_curve/s_curve_umap.rds")

scurve_scaled_obj <- gen_scaled_data(
  data = umap_scurve)

umap_scurve_scaled <- scurve_scaled_obj$scaled_nldr
lim1 <- scurve_scaled_obj$lim1
lim2 <- scurve_scaled_obj$lim2
r2 <- diff(lim2)/diff(lim1)

sc_ltr_pos <- c(0.08, 0.93)
# sc_xlims <- c(-0.5, 1.47)
# sc_ylims <- c(-0.32, 2.1)
sc_xlims <- c(-0.25, 1.25)
sc_ylims <- c(-0.3, 2)

nldr_scurve <- umap_scurve_scaled |>
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha=0.5, colour="#000000", size = 0.5) +
  xlim(sc_xlims) + ylim(sc_ylims) +
  interior_annotation("a", sc_ltr_pos)
```

<!--Full hexagon grid with UMAP data-->

  ```{r}
#| label: hexbin-scurve
## Compute hexbin parameters
num_bins_x_scurve <- 7

## hexagon binning to have regular hexagons
hb_obj_scurve <- hex_binning(
  data = umap_scurve_scaled,
  bin1 = num_bins_x_scurve,
  r2 = r2,
  q = 0.1)

## Data set with all centroids
all_centroids_df <- hb_obj_scurve$centroids

## Generate all coordinates of hexagons
hex_grid <- hb_obj_scurve$hex_poly

## To obtain the standardise counts within hexbins
counts_df <- hb_obj_scurve$std_cts
df_bin_centroids <- extract_hexbin_centroids(
  centroids_df = all_centroids_df,
  counts_df = counts_df) |>
  filter(drop_empty == FALSE)

hex_grid_with_counts <-
  left_join(hex_grid,
            counts_df,
            by = c("hex_poly_id" = "hb_id"))

hex_grid_scurve <- ggplot(
  data = hex_grid_with_counts,
  aes(x = x, y = y)) +
  geom_polygon(color = "black",
               aes(group = hex_poly_id),
               fill = "#ffffff") +
  geom_point(data = umap_scurve_scaled,
             aes(x = UMAP1, y = UMAP2),
             alpha = 0.5, size = 0.5, color = "#000000") +
  xlim(sc_xlims) + ylim(sc_ylims) +
  interior_annotation("b", sc_ltr_pos)
```

<!--Non-empty bins with bin centroids-->

  ```{r}
#| label: empty-bin-scurve
hex_grid_nonempty <- hex_grid |>
  filter(hex_poly_id %in% df_bin_centroids$hexID)

hex_grid_nonempty_scurve <- ggplot(
  data = hex_grid_nonempty,
  aes(x = x, y = y)) +
  geom_polygon(color = "black",
               aes(group = hex_poly_id),
               fill = "#ffffff") +
  geom_point(data = df_bin_centroids,
             aes(x = c_x, y = c_y),
             color = "#33a02c",
             size = 1) +
  xlim(sc_xlims) + ylim(sc_ylims) +
  interior_annotation("c", sc_ltr_pos)
```

<!--2D model-->

  ```{r}
#| label: triangulate-scurve
## Triangulate bin centroids
tr1_object_scurve <- tri_bin_centroids(
  df_bin_centroids, x = "c_x", y = "c_y")
tr_from_to_df_scurve <- gen_edges(
  tri_object = tr1_object_scurve)

## Compute 2D distances
distance_scurve <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_scurve,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_scurve <- find_lg_benchmark(
  distance_edges = distance_scurve,
  distance_col = "distance")

trimesh_removed_scurve <- vis_rmlg_mesh(
  distance_edges = distance_scurve,
  benchmark_value = benchmark_scurve,
  tr_coord_df = tr_from_to_df_scurve,
  distance_col = "distance") +
  xlim(sc_xlims) + ylim(sc_ylims) +
  interior_annotation("d", sc_ltr_pos)
```
```{r}
#| label: s-curve-true-model-proj1

true_model_df <- read_rds("data/s_curve/scurve_true_model.rds")
wireframe_true_model <- read_rds("data/s_curve/scurve_true_model_wireframe.rds")

data_scurve <- training_data_scurve |>
  select(-ID) |>
  mutate(type = "data")

model_scurve <- true_model_df |>
  select(-ID) |>
  mutate(type = "model")

# Apply the scaling
df_model_data_scurve <- bind_rows(data_scurve, model_scurve)
scaled_scurve <- scale_data_manual(df_model_data_scurve, "type") |>
  as_tibble()

scaled_s_curve_noise <- scaled_scurve |>
  filter(type == "data") |>
  select(-type)

scaled_s_curve_true <- scaled_scurve |>
  filter(type == "model") |>
  select(-type)

projection <- cbind(
  c(0.45181,0.04747,0.02613,0.25735,0.20622,0.08755,-0.11488),
  c(0.03231,0.23556,0.35166,-0.11741,0.12730,0.09544,0.34262))

projected <- as.matrix(scaled_s_curve_noise) %*% projection

projected_df <- projected |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::rename(c("proj1" = "...1",
                  "proj2" = "...2")) |>
  dplyr::mutate(ID = dplyr::row_number())

projected_true <- as.matrix(scaled_s_curve_true) %*% projection

projected_true_df <- projected_true |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::rename(c("proj1" = "...1",
                  "proj2" = "...2")) |>
  dplyr::mutate(ID = dplyr::row_number())


model_df <- dplyr::left_join(
  wireframe_true_model,
  projected_true_df,
  by = c("from" = "ID"))

names(model_df)[3:NCOL(model_df)] <- paste0(names(projected_true_df)[-NCOL(projected_true_df)], "_from")

model_df <- dplyr::left_join(model_df, projected_true_df, by = c("to" = "ID"))
names(model_df)[(2 + NCOL(projected_true_df)):NCOL(model_df)] <- paste0(names(projected_true_df)[-NCOL(projected_true_df)], "_to")

limits <- 1
rng <- range(projected)
projected <- projected/max(abs(rng))
colnames(projected) <- c("P1", "P2")
projected <- data.frame(projected)
obs_labels <- as.character(1:nrow(scaled_s_curve_noise))

axis_scale <- limits/6
axis_pos <- -2/3 * limits

# axis_scale <- limits/4
# axis_pos <- -0.5 * limits

adj <- function(x) axis_pos + x * axis_scale
axes <- data.frame(x1 = adj(0),
                   y1 = adj(0),
                   x2 = adj(projection[, 1]),
                   y2 = adj(projection[, 2]))

axis_labels <- colnames(scaled_s_curve_noise)
rownames(axes) <- axis_labels

# axis_scale <- 0.15 * limits
# axis_pos <- -0.6 * limits
# adj <- function(x) axis_pos + x * axis_scale
theta <- seq(0, 2 * pi, length = 50)
circle <- data.frame(c1 = adj(cos(theta)),
                     c2 = adj(sin(theta)))

scurve_proj1_true_model <- projected_df |>
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
    color = "#3182bd") +
  geom_point(
    size = 0.8,
    alpha = 0.5,
    color = "#000000") +
  geom_segment(
    data=axes,
    aes(x=x1, y=y1, xend=x2, yend=y2),
    colour="grey70") +
  geom_text(
    data=axes,
    aes(x=x2, y=y2, label=rownames(axes)),
    colour="grey50",
    size = 3) +
  geom_path(
    data=circle,
    aes(x=c1, y=c2), colour="grey70") +
  coord_fixed()

```

```{r}
#| label: s-curve-true-model-proj2

projection <- cbind(
  c(0.03627,0.15991,0.39768,-0.21032,-0.14886,-0.17589,0.23305),
  c(-0.34731,-0.08183,0.15819,0.31926,0.08105,0.10083,0.25627))

projected <- as.matrix(scaled_s_curve_noise) %*% projection

projected_df <- projected |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::rename(c("proj1" = "...1",
                  "proj2" = "...2")) |>
  dplyr::mutate(ID = dplyr::row_number())

projected_true <- as.matrix(scaled_s_curve_true) %*% projection

projected_true_df <- projected_true |>
  tibble::as_tibble(.name_repair = "unique") |>
  dplyr::rename(c("proj1" = "...1",
                  "proj2" = "...2")) |>
  dplyr::mutate(ID = dplyr::row_number())


model_df <- dplyr::left_join(
  wireframe_true_model,
  projected_true_df,
  by = c("from" = "ID"))

names(model_df)[3:NCOL(model_df)] <- paste0(names(projected_true_df)[-NCOL(projected_true_df)], "_from")

model_df <- dplyr::left_join(model_df, projected_true_df, by = c("to" = "ID"))
names(model_df)[(2 + NCOL(projected_true_df)):NCOL(model_df)] <- paste0(names(projected_true_df)[-NCOL(projected_true_df)], "_to")

limits <- 1
rng <- range(projected)
projected <- projected/max(abs(rng))
colnames(projected) <- c("P1", "P2")
projected <- data.frame(projected)
obs_labels <- as.character(1:nrow(scaled_s_curve_noise))

# axis_scale <- 2 * limits/3
# axis_pos <- -0.6 * limits
axis_scale <- limits/6
axis_pos <- -2/3 * limits

adj <- function(x) axis_pos + x * axis_scale
axes <- data.frame(x1 = adj(0),
                   y1 = adj(0),
                   x2 = adj(projection[, 1]),
                   y2 = adj(projection[, 2]))

axis_labels <- colnames(scaled_s_curve_noise)
rownames(axes) <- axis_labels

# axis_scale <- 0.15 * limits
# axis_pos <- -0.6 * limits
# adj <- function(x) axis_pos + x * axis_scale
theta <- seq(0, 2 * pi, length = 50)
circle <- data.frame(c1 = adj(cos(theta)), c2 = adj(sin(theta)))

scurve_proj2_true_model <- projected_df |>
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
    color = "#3182bd") +
  geom_point(
    size = 0.8,
    alpha = 0.5,
    color = "#000000") +
  geom_segment(
    data=axes,
    aes(x=x1, y=y1, xend=x2, yend=y2),
    colour="grey70") +
  geom_text(
    data=axes,
    aes(x=x2, y=y2, label=rownames(axes)),
    colour="grey50",
    size = 3) +
  geom_path(
    data=circle,
    aes(x=c1, y=c2), colour="grey70") +
  coord_fixed()

```

```{r}
#| label: fig-scurve-true-sc
#| fig-cap: "Two views of the true model (blue lines) in $2\\text{-}D$ projections from $7\\text{-}D$, for the S-curve data (black points). The data is spread along the S shape, and does not vary much from this curve. (The **langevitour** software is used to view the data with a tour, and the full video is available at <https://youtu.be/I5GL23vLiw0>)."
#| fig-pos: H
#| fig-width: 10
#| fig-height: 5

scurve_proj1_true_model + scurve_proj2_true_model +
  plot_layout(ncol = 2)
```

```{r}
#| label: fig-NLDR-scurve
#| echo: false
#| fig-cap: "Key steps for constructing the model on the UMAP layout ($k=2$): (a) data, (b) hexagon bins, (c) bin centroids, and (d) triangulated centroids. The S-curve data is shown."
#| fig-width: 8
#| fig-height: 4
#| out-width: 100%

nldr_scurve + hex_grid_scurve +
  hex_grid_nonempty_scurve +
  trimesh_removed_scurve +
  plot_layout(ncol = 4)
```

```{r}
#| label: code-illustration
# Code to draw illustration for notation
start_pt <- all_centroids_df |>
  filter(hexID == 1)
d_rect <- tibble(x1min = 0,
                 x1max = 1,
                 x2min = 0,
                 x2max = r2)

# To move the rectangle to ignore the overlap with the centroids
rect_adj <- tibble(x1 = 0.03, x2 = 0.03)

a1 <- tibble(x = all_centroids_df$c_x[4],
             xend = all_centroids_df$c_x[5],
             y = all_centroids_df$c_y[21],
             yend = all_centroids_df$c_y[21],
             label = expression(a[1]))
a2 <- tibble(x = all_centroids_df$c_x[25],
             xend = all_centroids_df$c_x[25],
             y = all_centroids_df$c_y[25],
             yend = all_centroids_df$c_y[33],
             label = expression(a[2]))
hex_param_vis <- ggplot() +
  geom_polygon(data = hex_grid,
               aes(x = x,
                   y = y,
                   group = hex_poly_id),
               fill = "white",
               color = "#bdbdbd") +
  geom_point(data = all_centroids_df, aes(
    x = c_x,
    y = c_y),
    color = "#31a354", size = 0.9) +
  geom_point(data = start_pt, aes(x = c_x,
                                  y = c_y),
             color = "black") +
  geom_rect(data=d_rect,
            aes(xmin = x1min - rect_adj$x1,# - rect_adj$s1,
                xmax = x1max - rect_adj$x1,# - rect_adj$s1,
                ymin = x2min - rect_adj$x2,# - rect_adj$s2,
                ymax = x2max - rect_adj$x2),# - rect_adj$s2),
            fill = "white",
            color = "black",
            alpha = 0,
            linewidth = 0.7) +
  geom_point(data=d_rect, aes(x=x1min - rect_adj$x1,
                              y=x2min - rect_adj$x2)) +
  geom_point(data=d_rect, aes(x=x1max - rect_adj$x1,
                              y=x2min - rect_adj$x2)) +
  geom_point(data=d_rect, aes(x=x1min - rect_adj$x1,
                              y=x2max - rect_adj$x2)) +
  annotate("text", x=d_rect$x1min - rect_adj$x1,
           y=d_rect$x2min - rect_adj$x2,
           label = "(0,0)",
           hjust=-0.1, vjust=-0.3) +
  annotate("text", x=d_rect$x1max - rect_adj$x1,
           y=d_rect$x2min - rect_adj$x2,
           label = "(0,1)",
           hjust=1.1, vjust=-0.3) +
  annotate("text", x=d_rect$x1min - rect_adj$x1,
           y=d_rect$x2max - rect_adj$x2,
           label = expression(group("(",
                                    list(0, y[2][max]),")")),
           hjust=-0.1, vjust=1.2) +
  geom_segment(data=d_rect, aes(
    x = x1min  - rect_adj$x1, # 0 - 0.03,
    y = -0.35,
    xend = x1max - rect_adj$x1, #1 - 0.03,
    yend = -0.35), #-0.35),
    arrow = arrow(length = unit(0.03, "npc"),
                  ends = "both"),
    color = "black")+
  annotate("text", x=0.5, y=-0.45,
           label = expression(r[1]), color = "black") +
  geom_segment(data=d_rect, aes(
    x = -0.25,
    y = x2min - rect_adj$x2, #0 - 0.05,
    xend = -0.25,
    yend = x2max - rect_adj$x2), #r2 - 0.05),
    arrow = arrow(length = unit(0.03, "npc"),
                  ends = "both"),
    color = "black")+
  annotate("text", x=-0.35, y=1,
           label = expression(r[2]), color = "black") +
  geom_segment(data = a1, aes(
    x = x, #-0.1 + 0.2087578,
    y = y, #-0.15,
    xend = xend, #-0.1 + 0.2087578*2,
    yend = yend), #-0.15),
    arrow = arrow(length = unit(0.03, "npc"),
                  ends = "both"),
    color = "black")+ # a1 = 0.2087578
  annotate("text",
           x=(a1$x+a1$xend)/2,
           y=a1$y,
           label = expression(a[1]),
           color = "black",
           vjust = 1.2) +
  geom_segment(data = a2, aes(
    x = x, #-0.15,
    y = y, #-0.1*r2 + 0.1807896*2,
    xend = xend, #-0.15,
    yend = yend), #-0.1*r2 + 0.1807896*3),
    arrow = arrow(length = unit(0.03, "npc"),
                  ends = "both"),
    color = "black") + # a2 = 0.1807896
  annotate("text", x=a2$x, y=(a2$y+a2$yend)/2,
           label = expression(a[2]),
           color = "black", hjust=-0.2) +
  annotate("text", x=-0.18, y=-0.25,
           label = expression(group("(", list(s[1], s[2]), ")")),
           color = "black") +
  coord_equal()
```

```{r}
#| label: fig-hex-param
#| fig-cap: "The components of the hexagon grid illustrating notation."
#| out-height: 30%
#| fig-pos: H

hex_param_vis
```

```{r}
#| label: best-umap-model-scurve

# Initialize effective bins along x
effective_bin1_scurve <- 15

# Fit high-dimensional model
scurve_model <- fit_highd_model(
  training_data = training_data_scurve,
  emb_df = umap_scurve_scaled,
  bin1 = effective_bin1_scurve,
  r2 = r2,
  q = 0.1,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_scurve <- scurve_model$df_bin_centroids
df_bin_scurve <- scurve_model$df_bin

hex_grid_nonempty <- hex_grid |>
  filter(hex_poly_id %in% df_bin_centroids_scurve$hexID)

```

```{r}
#| label: p-d-error-in-2d-scurve

## Compute error
error_df_scurve <- augment(
  df_bin_centroids = df_bin_centroids_scurve,
  df_bin = df_bin_scurve,
  training_data = training_data_scurve,
  newdata = NULL,
  type_NLDR = "UMAP",
  col_start = "x")

# error_df_scurve <- error_df_scurve |>
#   mutate(type = case_when(
#     row_wise_abs_error == 0 ~ "no error",
#     row_wise_abs_error <= 0.2 ~ "error 0-0.2",
#     row_wise_abs_error <= 0.4 ~ "error 0.2-0.4",
#     row_wise_abs_error <= 0.6 ~ "error 0.4-0.6",
#     .default = "error greter than 0.6"
#   )) |>
#   mutate(type = factor(type, levels = c(
#     "no error", "error 0-0.2", "error 0.2-0.4", "error 0.4-0.6", "error greter than 0.6")))


error_df_scurve <- error_df_scurve |>
  mutate(sqrt_row_wise_abs_error = sqrt(row_wise_abs_error))

error_df_scurve <- error_df_scurve |>
  bind_cols(umap_scurve_scaled |>
              select(-ID))

error_plot_scurve <- error_df_scurve |>
  ggplot(aes(x = UMAP1,
             y = UMAP2,
             colour = sqrt_row_wise_abs_error)) +
  geom_point(alpha=0.5) +
  scale_colour_continuous_sequential(palette = "Light Grays") +
  coord_equal()
```

```{r}
#| label: fig-p-d-error-in-2d-scurve
#| fig-cap: "error"
#| out-height: 30%
#| fig-pos: H

error_plot_scurve
```

```{r}
#| label: model-prediction-umap-original

predict_umap_scurve <- read_rds(file = "data/s_curve/s_curve_umap_predict_true.rds")

predict_scurve_obj <- gen_scaled_data(
  data = predict_umap_scurve)

predict_umap_scurve_scaled <- predict_scurve_obj$scaled_nldr

plot_predict_umap <- ggplot(data = predict_umap_scurve_scaled,
                            aes(
                              x = UMAP1,
                              y = UMAP2
                            )) +
  geom_point(alpha = 0.5) +
  interior_annotation("b", c(0.92, 0.96))
```

```{r}
#| label: model-prediction-umap-quollr

true_pred_df <- predict_emb(
  test_data = true_model_df,
  df_bin_centroids = df_bin_centroids_scurve,
  df_bin = df_bin_scurve,
  type_NLDR = "UMAP"
)


# Compute radius r
a1 <- calc_bins_y(bin1 = 15, r2 = r2, q = 0.16)$a1
r <- a1 / 2

# Function to jitter points within a circumcircle of a hexagon
jitter_within_circumcircle <- function(center_x, center_y, radius, num_points) {
  theta <- runif(num_points, 0, 2 * pi)
  r <- radius * sqrt(runif(num_points))
  x <- center_x + r * cos(theta)
  y <- center_y + r * sin(theta)
  jittered_points <- data.frame(UMAP1 = x, UMAP2 = y, ID = points_in_hex$ID)
  return(jittered_points)
}

# Jittered points data frame
jittered_points_df <- data.frame()

# Iterate through each hexagon and jitter points within the circumcircle
for (hex_id in unique(true_pred_df$pred_hb_id)) {
  # hex_points <- hex_grid_nonempty %>% filter(hex_poly_id == hex_id)
  # center_x <- mean(hex_points$x)
  # center_y <- mean(hex_points$y)
  # points_in_hex <- true_pred_df %>% filter(pred_UMAP_1 >= min(hex_points$x) & pred_UMAP_1 <= max(hex_points$x) & pred_UMAP_2 >= min(hex_points$y) & pred_UMAP_2 <= max(hex_points$y))

  points_in_hex <- true_pred_df |>
    filter(pred_hb_id == hex_id)

  center_x <- points_in_hex |>
    pull(pred_UMAP_1) |>
    unique()

  center_y <- points_in_hex |>
    pull(pred_UMAP_2) |>
    unique()

  if (nrow(points_in_hex) > 0) {
    jittered_points <- jitter_within_circumcircle(center_x, center_y, r, nrow(points_in_hex))
    jittered_points_df <- rbind(jittered_points_df, jittered_points)
  }
}

plot_predict_umap_model <- ggplot(data = jittered_points_df,
                                  aes(
                                    x = UMAP1,
                                    y = UMAP2
                                  )) +
  geom_point(alpha = 0.5) +
  interior_annotation("c", c(0.92, 0.96))
```

```{r}
#| label: true-model-projection-langevitour

true_model_df_rm_id <- true_model_df |>
  select(-ID)

projection <- cbind(
  c(0.18390,-0.18443,-0.06860,0.00654,-0.02921,0.02313,-0.08963),
  c(-0.02902,0.03638,-0.16872,0.22418,0.01146,0.02777,0.01453))
projected_df <- as.matrix(true_model_df_rm_id) %*% projection |>
  as_tibble(.name_repair = "unique") |>
  rename(c("proj1" = "...1",
           "proj2" = "...2"))

proj_plot_scurve_true <- ggplot(
  data = projected_df,
  aes(
    x = proj1,
    y = proj2)) +
  geom_point(alpha = 0.5) +
  interior_annotation("a", c(0.92, 0.96))
```

```{r}
#| label: model-prediction-points-line

predict_umap_scurve_scaled <- predict_umap_scurve_scaled |>
  mutate(type = "predict_from_umap")

jittered_points_df <- jittered_points_df |>
  mutate(type = "predict_from_model")

predict_point_df <- bind_rows(
  predict_umap_scurve_scaled,
  jittered_points_df
)

plot_predict_umap_model_connect <- ggplot(data = predict_point_df,
                                          aes(
                                            x = UMAP1,
                                            y = UMAP2,
                                            group = ID
                                          )) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.3) +
  interior_annotation("d", c(0.92, 0.96))
```

```{r}
#| echo: false
#| label: fig-predict-scurve
#| fig-pos: H
#| fig-width: 9
#| fig-height: 5
#| fig-cap: "Comparison of prediction generated using the exiting `umap` prediction method and our method: (a) A view of the true model in projections from $7\\text{-}D$, (b) predicted data from the `umap` prediction method, (c) predicted data from our method, and (d) comparison of predictions from both methods. In plot (d), points representing the same data but predicted by different methods are connected by lines. Need to add concluded sentence after changing the data."

proj_plot_scurve_true + plot_predict_umap +
  plot_predict_umap_model + plot_predict_umap_model_connect +
  plot_layout(ncol = 4)
```

```{r}
#| label: hexbin-regular-scurve1
## hexagon binning to have regular hexagons
hb_obj_scurve1 <- hex_binning(
  data = umap_scurve_scaled,
  bin1 = 7,
  r2 = r2)

a1_1 <- calc_bins_y(
  bin1 = 7,
  r2 = r2
)$a1

## Data set with all centroids
all_centroids_df1 <- hb_obj_scurve1$centroids

## Generate all coordinates of hexagons
hex_grid1 <- hb_obj_scurve1$hex_poly

## To obtain the standardise counts within hexbins
counts_df1 <- hb_obj_scurve1$std_cts
df_bin_centroids1 <- extract_hexbin_centroids(
  centroids_df = all_centroids_df1,
  counts_df = counts_df1) |>
  filter(drop_empty == FALSE) |>
  mutate(b1 = "b1 = 7")

hex_grid_with_counts_s_curve1 <- full_join(
  hex_grid1,
  df_bin_centroids1 |> select(hexID, std_counts),
  by = c("hex_poly_id" = "hexID"))

hex_grid_coloured_scurve1 <- ggplot() +
  geom_polygon(
    data = hex_grid_with_counts_s_curve1,
    aes(x = x, y = y,
        group = hex_poly_id,
        fill = std_counts), color = "grey70", linewidth=0.2) +
  geom_point(data = umap_scurve_scaled,
             aes(x = UMAP1, y = UMAP2),
             alpha = 0.3,
             size = 0.5) +
  scale_fill_viridis_c(direction = -1,
                       na.value = "#ffffff", option = "E") +
  interior_annotation("a", position = c(0.92, 0.93)) +
  xlim(c(-0.25, 1.35)) + ylim(c(-0.35, 2))

hex_grid_coloured_scurve1_dens <- ggplot(hex_grid_with_counts_s_curve1) +
  geom_histogram(aes(x=std_counts),
                 breaks=seq(0, 1, 0.05),
                 fill="slategray4", colour="white")

#hex_grid_coloured_scurve1_dens <-
#  hex_grid_with_counts_s_curve1 |>
#  filter(!is.na(std_counts)) |>
#  mutate(ord_std_counts = n() - rank(std_counts)) |>
#  ggplot() +
#  geom_line(aes(x=ord_std_counts, y=std_counts),
#                 colour="slategray4")
```

```{r}
#| label: hexbin-regular-scurve2
## hexagon binning to have regular hexagons
hb_obj_scurve2 <- hex_binning(
  data = umap_scurve_scaled,
  bin1 = 11,
  r2 = r2)

a1_2 <- calc_bins_y(
  bin1 = 11,
  r2 = r2
)$a1

## Data set with all centroids
all_centroids_df2 <- hb_obj_scurve2$centroids

## Generate all coordinates of hexagons
hex_grid2 <- hb_obj_scurve2$hex_poly

## To obtain the standardise counts within hexbins
counts_df2 <- hb_obj_scurve2$std_cts
df_bin_centroids2 <- extract_hexbin_centroids(
  centroids_df = all_centroids_df2,
  counts_df = counts_df2) |>
  filter(drop_empty == FALSE) |>
  mutate(b1 = "b1 = 11")

hex_grid_with_counts_s_curve2 <- full_join(
  hex_grid2,
  df_bin_centroids2 |> select(hexID, std_counts),
  by = c("hex_poly_id" = "hexID"))

hex_grid_coloured_scurve2 <- ggplot() +
  geom_polygon(
    data = hex_grid_with_counts_s_curve2,
    aes(x = x, y = y,
        group = hex_poly_id,
        fill = std_counts),
    color = "grey70",
    linewidth=0.2) +
  geom_point(data = umap_scurve_scaled,
             aes(x = UMAP1, y = UMAP2),
             alpha = 0.3,
             size = 0.5) +
  scale_fill_viridis_c(direction = -1,
                       na.value = "#ffffff", option = "C") +
  xlim(c(-0.25, 1.35)) + ylim(c(-0.35, 2)) +
  interior_annotation("b")

hex_grid_coloured_scurve2_dens <- ggplot(hex_grid_with_counts_s_curve2) +
  geom_histogram(aes(x=std_counts),
                 breaks=seq(0, 1, 0.05),
                 fill="slategray4", colour="white")

#hex_grid_coloured_scurve2_dens <-
#  hex_grid_with_counts_s_curve2 |>
#  filter(!is.na(std_counts)) |>
#  mutate(ord_std_counts = n() - rank(std_counts)) |>
#  ggplot() +
#  geom_line(aes(x=ord_std_counts, y=std_counts),
#                 colour="slategray4")

```

```{r}
#| label: hexbin-regular-scurve3
## hexagon binning to have regular hexagons
hb_obj_scurve3 <- hex_binning(
  data = umap_scurve_scaled,
  bin1 = 13,
  r2 = r2)

a1_3 <- calc_bins_y(
  bin1 = 13,
  r2 = r2
)$a1

## Data set with all centroids
all_centroids_df3 <- hb_obj_scurve3$centroids

## Generate all coordinates of hexagons
hex_grid3 <- hb_obj_scurve3$hex_poly

## To obtain the standardise counts within hexbins
counts_df3 <- hb_obj_scurve3$std_cts
df_bin_centroids3 <- extract_hexbin_centroids(
  centroids_df = all_centroids_df3,
  counts_df = counts_df3) |>
  filter(drop_empty == FALSE) |>
  mutate(b1 = "b1 = 14")

hex_grid_with_counts_s_curve3 <- full_join(hex_grid3,
                                           df_bin_centroids3 |> select(hexID, std_counts),
                                           by = c("hex_poly_id" = "hexID"))

hex_grid_coloured_scurve3 <-  ggplot() +
  geom_polygon(
    data = hex_grid_with_counts_s_curve3,
    aes(x = x, y = y,
        group = hex_poly_id,
        fill = std_counts),
    color = "grey70", linewidth=0.2) +
  geom_point(data = umap_scurve_scaled,
             aes(x = UMAP1, y = UMAP2),
             alpha = 0.3,
             size = 0.5) +
  scale_fill_viridis_c(direction = -1, na.value = "#ffffff", option = "D") +
  xlim(c(-0.25, 1.35)) + ylim(c(-0.35, 2)) +
  interior_annotation("c")

hex_grid_coloured_scurve3_dens <- ggplot(hex_grid_with_counts_s_curve3) +
  geom_histogram(aes(x=std_counts),
                 breaks=seq(0, 1, 0.05),
                 fill="slategray4", colour="white")
#hex_grid_coloured_scurve3_dens <-
#  hex_grid_with_counts_s_curve3 |>
#  filter(!is.na(std_counts)) |>
#  mutate(ord_std_counts = n() - rank(std_counts)) |>
#  ggplot() +
#  geom_line(aes(x=ord_std_counts, y=std_counts),
#                 colour="slategray4")
```

```{r}
#| echo: false
#| label: fig-bins-scurve
#| fig-pos: H
#| fig-cap: "Hexbin density plots of UMAP layout of the S-curve data, using three different bin inputs: (a) $b = 91 \\text{ } (7, \\text{ }13)$, (b) $b = 220 \\text{ } (11, \\text{ }20)$, and (c) $b = 312 \\text{ } (13, \\text{ }24)$. Color indicates standardized counts, dark indicating high count and light indicates low count. At the smallest bin size the data segregates into two separate groups, suggesting this is too many bins. Using the MSE of the model fit in $p\\text{-}D$ helps decide on a useful choice of number of bins."
#| fig-width: 6
#| fig-height: 4

hex_grid_coloured_scurve1 +
  hex_grid_coloured_scurve2 +
  hex_grid_coloured_scurve3 +
  plot_layout(guides='collect', ncol = 3,
              heights = c(2,1)) &
  theme(legend.position='none', plot.tag = element_text(size = 8))
```

The default number of bins $b=b_1\times b_2$ is computed based on the sample size, by setting $b_1=n^{1/3}$, consistent with the Diaconis-Freedman rule [@freedman1981]. The value of $b_2$ is determined analytically by $b_1, q, r_2$. Values of $b_1$ between $2$ and $b_1 = \sqrt{\frac{n}{r_2}}$ are allowed. @fig-param-scurve (b) shows the effect of different choices of $b_1$ on the MSE of the fitted model.

<!-- To determine the effective $b$, candidate values are selected based on the range between the minimum and approximate maximum $b_1$, because $b_2$ is computed from $b_1$. The minimum $b_1$ is set to $2$, while the maximum number is estimated by taking the square root of $\frac{n}{2}$. By evaluating MSE across varying $b$ within this range for different $q$, helps to determine an appropriate values for $b$ and $q$ (@fig-param-scurve (a)).-->

  <!--add MSE vs total number of error plot-->

  <!--To generate errors for different total number of bins-->

  ```{r}
#| label: errors-scurve
## To initialize number of bins along the x-axis
bin1_vec_scurve <- 5:22 #sqrt(NROW(umap_scurve_scaled)/r2)
#buffer_vec <- seq(0.05, 0.1, by = 0.01)

error_scurve <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_scurve) {
  #for (buff in buffer_vec) {

  bin_obj <- calc_bins_y(
    bin1 = xbins,
    r2 = r2,
    q = 0.1)

  bin2 <- bin_obj$bin2
  a1 <- bin_obj$a1
  a2 <- bin_obj$a2

  scurve_model <- fit_highd_model(
    training_data = training_data_scurve,
    emb_df = umap_scurve_scaled,
    bin1 = xbins,
    r2 = r2,
    q = 0.1,
    is_bin_centroid = TRUE,
    is_rm_lwd_hex = FALSE,
    col_start_highd = "x"
  )

  df_bin_centroids_scurve <- scurve_model$df_bin_centroids
  df_bin_scurve <- scurve_model$df_bin

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_scurve,
    df_bin = df_bin_scurve,
    training_data = training_data_scurve,
    newdata = NULL,
    type_NLDR = "UMAP",
    col_start = "x") |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_scurve),
           a1 = round(a1, 2),
           a2 = round(a2, 2),
           side_length = quad(a=3, b = 2 * a2, c = -(a2^2 + a1^2)))

  error_scurve <- bind_rows(error_scurve, error_df)

  #}

}

# mse_scurve_b <- ggplot(error_scurve,
#                      aes(x = b,
#                          y = log(MSE),
#                          color = buff,
#                          group = buff)) +
#   geom_point(size = 1) +
#   geom_line(linewidth = 0.05) +
#   geom_vline(xintercept = 220, linetype="solid",
#              color = "black", linewidth=0.1) +
#   geom_vline(xintercept = 112, linetype=2,
#              color = "black", linewidth=0.1) +
#   geom_vline(xintercept = 364, linetype=2,
#              color = "black", linewidth=0.1) +
#   scale_colour_continuous_sequential(palette = "Magenta") +
#   labs(x = expression(b), y = "log(MSE)", color = expression(q)) +
#   ggtitle("(a)") +
#   theme(aspect.ratio = 0.75,
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12,
#             angle=90),
#         plot.title = element_text(size = 12),
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.key.height = unit(1, 'cm'),
#         legend.key.width = unit(1, 'cm'),
#         legend.key.size = unit(1, 'cm'))

# diaconis_bin_width <- calc_bins_y(
#     bin1 = ceiling(NROW(training_data_scurve)^(1/3)),
#     r2 = r2,
#     q = 0.1)$a1


## Find the minimum MSE when have duplicate a1
error_scurve <- error_scurve |>
  group_by(a1) |>
  filter(MSE == min(MSE)) |>
  ungroup() |>
  mutate(prop_dens = b_non_empty/side_length)

base_line_prop_dens <- error_scurve |>
  filter(a1 == min(a1)) |>
  pull(prop_dens)

error_scurve <- error_scurve |>
  mutate(prop_comp = prop_dens/base_line_prop_dens)

mse_scurve_b <- ggplot(error_scurve,
                       aes(x = a1,
                           y = MSE)) +
  geom_vline(xintercept = 0.10, linetype=2,
             color = "#54278f", linewidth=1) +
  geom_vline(xintercept = 0.12,
             linetype="solid",
             color = "#006d2c",
             linewidth=1) +
  geom_vline(xintercept = 0.18, linetype=2,
             color = "#08519c", linewidth=1) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  scale_x_continuous(breaks = sort(unique(round(error_scurve$a1, 2)))[c(1, 3, 5, 6, 8, 10:12)]) +
  labs(x = expression(paste("binwidth (", a[1], ")")), y = "MSE") +
  ggtitle("(a)") +
  theme_minimal() +
  theme(aspect.ratio = 0.75,
        panel.border = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 12, hjust = 0.5, vjust = -0.5),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_line())

prop_dens_a1 <- ggplot(error_scurve,
                       aes(x = a1,
                           y = prop_comp)) +
  geom_vline(xintercept = 0.12, linetype="solid",
             color = "#006d2c", linewidth=1) +
  geom_vline(xintercept = 0.10, linetype=2,
             color = "#54278f", linewidth=1) +
  geom_vline(xintercept = 0.18, linetype=2,
             color = "#08519c", linewidth=1) +
  geom_point(size = 1) +
  geom_line(linewidth = 0.05) +
  # scale_y_continuous(breaks = sort(unique(error_scurve$b_non_empty))[-3]) +
  scale_x_continuous(breaks = sort(unique(round(error_scurve$a1, 2)))[c(1, 3, 5, 6, 8, 10:12)]) +
  labs(x = expression(paste("binwidth (", a[1], ")")), y = paste("proportion of densities")) +
  ggtitle("(d)") +
  theme_minimal() +
  theme(aspect.ratio = 0.75,
        panel.border = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 12, hjust = 0.5, vjust = -0.5),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_line())


error_scurve <- error_scurve |>
  mutate(prop_bins = b_non_empty/b)

a1_m_scurve <- ggplot(error_scurve,
                      aes(x = a1,
                          y = prop_bins)) +
  geom_vline(xintercept = 0.12, linetype="solid",
             color = "#006d2c", linewidth=1) +
  geom_vline(xintercept = 0.10, linetype=2,
             color = "#54278f", linewidth=1) +
  geom_vline(xintercept = 0.18, linetype=2,
             color = "#08519c", linewidth=1) +
  geom_point(size = 1) +
  geom_line(linewidth = 0.05) +
  # scale_y_continuous(breaks = sort(unique(error_scurve$b_non_empty))[-3]) +
  scale_x_continuous(breaks = sort(unique(round(error_scurve$a1, 2)))[c(1, 3, 5, 6, 8, 10:12)]) +
  labs(x = expression(paste("binwidth (", a[1], ")")), y = expression(paste("proportion of non-empty bins ", bgroup("(", frac(m, b), ")")))) +
  ggtitle("(b)") +
  theme_minimal() +
  theme(aspect.ratio = 0.75,
        panel.border = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 12, hjust = 0.5, vjust = -0.5),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_line())

#
# b1_m_scurve2 <- ggplot(error_scurve,
#                      aes(x = bin1,
#                          y = b_non_empty)) +
#   geom_point(size = 1) +
#   geom_line(linewidth = 0.05) +
#   # geom_vline(xintercept = 11, linetype="solid",
#   #            color = "#bdbdbd", linewidth=1) +
#   # geom_vline(xintercept = 9, linetype=2,
#   #            color = "#bdbdbd", linewidth=1) +
#   # geom_vline(xintercept = 13, linetype=2,
#   #            color = "#bdbdbd", linewidth=1) +
#   scale_x_continuous(breaks = 5:19) +
#   labs(x = expression(b[1]), y = expression(m)) +
#   ggtitle("(c)") +
#   theme(aspect.ratio = 0.75,
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12,
#             angle=90),
#         axis.ticks.x = element_line(linewidth = 0.5),
#         axis.ticks.y = element_line(linewidth = 0.5),
#         plot.title = element_text(size = 12),
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.key.height = unit(1, 'cm'),
#         legend.key.width = unit(1, 'cm'),
#         legend.key.size = unit(1, 'cm'))
```

```{r}
#| label: distribution-scurve
## To initialize effective bins along x
effective_bin1_scurve <- 11

effective_bin2_scurve <- calc_bins_y(
  bin1 = effective_bin1_scurve,
  r2 = r2)$bin2

scurve_model <- fit_highd_model(
  training_data = training_data_scurve,
  emb_df = umap_scurve_scaled,
  bin1 = effective_bin1_scurve,
  r2 = r2,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_scurve <- scurve_model$df_bin_centroids
df_bin_scurve <- scurve_model$df_bin

## To obtain the first quntile
benchmark1 <- round(quantile(df_bin_centroids_scurve$std_counts, names = FALSE)[2], 3)
```

```{r}
#| label: bin-counts-scurve
#| eval: false

df_bin_centroids_all <- bind_rows(
  df_bin_centroids1,
  df_bin_centroids2,
  df_bin_centroids3
)

# cell_count_scurve <- ggplot(df_bin_centroids_all,
#                             aes(x = reorder(as.factor(hexID),
#                                             -bin_counts),
#                                 y = bin_counts,
#                                 group = as.factor(b1),
#                                 color = as.factor(b1))) +
#   geom_line() +
#   # geom_hline(yintercept = benchmark1, color = "black", linewidth=0.5, alpha = 0.5) +
#   scale_colour_discrete_qualitative(palette = "Dark 3") +
#   xlab("hexagon id (re-ordered)") +
#   ylab("bin count") +
#   ggtitle("(a)") +
#   theme(aspect.ratio = 0.75,
#         axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12, angle = 90),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 12),
#         plot.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.key.size = unit(1, 'cm'),
#         legend.key.height = unit(1, 'cm'))

cell_count_scurve <- ggplot(df_bin_centroids_all,
                            aes(x = as.factor(hexID),
                                y = std_counts,
                                group = as.factor(b1),
                                color = as.factor(b1))) +
  geom_line() +
  # geom_hline(yintercept = benchmark1, color = "black", linewidth=0.5, alpha = 0.5) +
  scale_colour_discrete_qualitative(palette = "Dark 3") +
  xlab("hexagon id (re-ordered)") +
  ylab("standardised bin count") +
  ggtitle("(a)") +
  theme(aspect.ratio = 0.75,
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'))
```

```{r}
#| label: bins-scurve

## To initialize benchmark values to remove low density hexagons
benchmark_rm_hex_vec <- sort(append(seq(0, 0.6, by=0.15), benchmark1))

## option 1
error_rm_scurve <- data.frame(matrix(nrow = 0, ncol = 0))

for (benchmark_rm_lwd in benchmark_rm_hex_vec) {

  df_bin_centroids_scurve_high_dens <- df_bin_centroids_scurve |>
    filter(std_counts > benchmark_rm_lwd)

  df_bin_scurve_high_dens <- df_bin_scurve |>
    filter(hb_id %in% df_bin_centroids_scurve_high_dens$hexID)

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_scurve_high_dens,
    df_bin = df_bin_scurve_high_dens,
    training_data = training_data_scurve,
    newdata = NULL,
    type_NLDR = "UMAP",
    col_start = "x") |>
    mutate(benchmark_rm_lwd = round(benchmark_rm_lwd, 2),
           bin1 = effective_bin1_scurve,
           bin2 = effective_bin2_scurve,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_scurve_high_dens),
           mean_counts = sum(df_bin_centroids_scurve_high_dens$bin_counts)/NROW(df_bin_centroids_scurve_high_dens))

  error_rm_scurve <- bind_rows(error_rm_scurve, error_df)

}

# benchmark_label_df <- tibble(x = 0.2, y = 0.1,
#                           label = paste0("benchmark is ", benchmark1))

error_rm_scurve <- error_rm_scurve |>
  mutate(mean_counts = round(mean_counts, 0)) |>
  filter(b_non_empty >= 5)

### Option 2

scurve_model2 <- fit_highd_model(
  training_data = training_data_scurve,
  emb_df = umap_scurve_scaled,
  bin1 = 7,
  r2 = r2,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_scurve2 <- scurve_model2$df_bin_centroids
df_bin_scurve2 <- scurve_model2$df_bin

error_rm_scurve2 <- data.frame(matrix(nrow = 0, ncol = 0))

for (benchmark_rm_lwd in benchmark_rm_hex_vec) {

  df_bin_centroids_scurve_high_dens <- df_bin_centroids_scurve2 |>
    filter(std_counts > benchmark_rm_lwd)

  df_bin_scurve_high_dens <- df_bin_scurve2 |>
    filter(hb_id %in% df_bin_centroids_scurve_high_dens$hexID)

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_scurve_high_dens,
    df_bin = df_bin_scurve_high_dens,
    training_data = training_data_scurve,
    newdata = NULL,
    type_NLDR = "UMAP",
    col_start = "x") |>
    mutate(benchmark_rm_lwd = round(benchmark_rm_lwd, 2),
           bin1 = 7,
           bin2 = 13,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_scurve_high_dens),
           mean_counts = sum(df_bin_centroids_scurve_high_dens$bin_counts)/NROW(df_bin_centroids_scurve_high_dens))

  error_rm_scurve2 <- bind_rows(error_rm_scurve2, error_df)

}

### Option 3

scurve_model3 <- fit_highd_model(
  training_data = training_data_scurve,
  emb_df = umap_scurve_scaled,
  bin1 = 13,
  r2 = r2,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_scurve3 <- scurve_model3$df_bin_centroids
df_bin_scurve3 <- scurve_model3$df_bin


error_rm_scurve3 <- data.frame(matrix(nrow = 0, ncol = 0))

for (benchmark_rm_lwd in benchmark_rm_hex_vec) {

  df_bin_centroids_scurve_high_dens <- df_bin_centroids_scurve3 |>
    filter(std_counts > benchmark_rm_lwd)

  df_bin_scurve_high_dens <- df_bin_scurve3 |>
    filter(hb_id %in% df_bin_centroids_scurve_high_dens$hexID)

  ## Compute error
  error_df <- glance(
    df_bin_centroids = df_bin_centroids_scurve_high_dens,
    df_bin = df_bin_scurve_high_dens,
    training_data = training_data_scurve,
    newdata = NULL,
    type_NLDR = "UMAP",
    col_start = "x") |>
    mutate(benchmark_rm_lwd = round(benchmark_rm_lwd, 2),
           bin1 = 13,
           bin2 = 24,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_scurve_high_dens),
           mean_counts = sum(df_bin_centroids_scurve_high_dens$bin_counts)/NROW(df_bin_centroids_scurve_high_dens))

  error_rm_scurve3 <- bind_rows(error_rm_scurve3, error_df)

}

error_rm_scurve <- bind_rows(
  error_rm_scurve,
  error_rm_scurve2,
  error_rm_scurve3
) |>
  mutate(bin1  = as.factor(bin1))

mse_scurve_lwd <- ggplot(error_rm_scurve,
                         aes(x = benchmark_rm_lwd,
                             y = MSE,
                             color = bin1)) +
  geom_point(
    size = 1
  ) +
  geom_line(
    linewidth = 0.3
  ) +
  # geom_vline(xintercept = benchmark1, linetype="solid",
  #            color = "#bdbdbd", linewidth=1, alpha = 0.5) +
  # scale_x_continuous("standardized bin count",
  #        transform = "reverse") +
  scale_x_continuous(breaks = unique(error_rm_scurve$benchmark_rm_lwd)) +
  scale_color_manual(values=c("#08519c", "#006d2c", "#54278f")) +
  xlab("threshold to remove low-density hexagons") +
  ylab("MSE") +
  ggtitle("(c)") +
  theme_minimal() +
  theme(aspect.ratio = 0.75,
        panel.border = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 12, hjust = 0.5, vjust = -0.5),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_line(),
        legend.position = "none")
```

```{r}
#| label: bin-neighbours
## First define low density bins using first quantile
df_bin_centroids_low <- df_bin_centroids_scurve |>
  filter(std_counts <= benchmark1)

## Check neighboring bins
remove_id <- find_low_dens_hex(df_bin_centroids_all = df_bin_centroids_scurve,
                               bin1 = effective_bin1_scurve,
                               df_bin_centroids_low = df_bin_centroids_low)

## Remove the identified bins
df_bin_centroids_scurve_removed <- df_bin_centroids_scurve |>
  filter(hexID %in%remove_id)

df_bin_centroids_scurve_keep  <- df_bin_centroids_scurve |>
  filter(!(hexID %in%remove_id))

```

```{r}
#| echo: false
#| fig-cap: "Various plots to help assess best hexagon bin parameters, thresholds to remove low density bins and large edges. Both (b) and (c) show MSE, against number of bins along the x-axis and standardised count. A good benchmark value for these parameters is when the MSE drops and then flattens out. Plot (a) shows the distribution of stadardised counts of hexagons. Plot (d) shows the distribution of $2\\text{-}D$ Euclidean distances between bin centroids, with a good benchmark value for removing large edges would being the distance that shows the first large decrease."
#| label: fig-param-scurve
#| out-width: 100%
#| fig-width: 12
#| fig-height: 10
#| fig-pos: H

mse_scurve_b + a1_m_scurve +
  mse_scurve_lwd + prop_dens_a1 +
  plot_layout(ncol=2)
```

<!--distribution of distance along with the default benchmark-->
  ```{r}
#| label: triangulate-scurve2

## Triangulate bin centroids
tr1_object_scurve <- tri_bin_centroids(
  df_bin_centroids_scurve_keep, x = "c_x", y = "c_y")
tr_from_to_df_scurve <- gen_edges(
  tri_object = tr1_object_scurve)

trimesh_scurve <- ggplot() +
  geom_segment(data = tr_from_to_df_scurve,
               aes(
                 x = x_from,
                 y = y_from,
                 xend = x_to,
                 yend = y_to),
               colour = "#33a02c") +
  interior_annotation("a", position = c(0.92, 0.95))

## Compute 2D distances
distance_scurve <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_scurve,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_scurve <- find_lg_benchmark(
  distance_edges = distance_scurve,
  distance_col = "distance")

## Benchmark 1

distance_df_small_edges_scurve1 <- distance_scurve |>
  filter(distance < benchmark_scurve)

tr_from_to_df_scurve1 <- inner_join(
  tr_from_to_df_scurve, distance_df_small_edges_scurve1,
  by = c("from", "to"))

trimesh_scurve_removed1 <- ggplot() +
  geom_segment(data = tr_from_to_df_scurve1,
               aes(
                 x = x_from,
                 y = y_from,
                 xend = x_to,
                 yend = y_to),
               colour = "#33a02c") +
  interior_annotation("c", position = c(0.92, 0.95))

## Benchmark 2

benchmark_scurve2 <- 0.6

distance_df_small_edges_scurve2 <- distance_scurve |>
  filter(distance < benchmark_scurve2)

tr_from_to_df_scurve2 <- inner_join(
  tr_from_to_df_scurve, distance_df_small_edges_scurve2,
  by = c("from", "to"))

trimesh_scurve_removed2 <- ggplot() +
  geom_segment(data = tr_from_to_df_scurve2,
               aes(
                 x = x_from,
                 y = y_from,
                 xend = x_to,
                 yend = y_to),
               colour = "#33a02c") +
  interior_annotation("b", position = c(0.92, 0.95))



## To draw the distance distribution
# distance_scurve$group <- "1"
# distance_scurve_plot <- ggplot(distance_scurve,
#                         aes(x = group,
#                             y = distance)) +
#   geom_quasirandom(
#     size = 2,
#     alpha = 0.3
#   ) +
#   ylim(0, max(unlist(distance_scurve$distance))+ 0.5) +
#   coord_flip() +
#   geom_hline(yintercept = benchmark_scurve,
#      linetype="solid",
#      color = "black", linewidth=0.5, alpha = 0.5) +
#   ylab(expression(d^{(2)})) +
#   ggtitle("(d)") +
#   theme(aspect.ratio = 0.75,
#         axis.text.x = element_text(size = 12),
#         axis.title.x = element_text(size = 12),
#         plot.title = element_text(size = 12))

# Define bin width and starting point
bin_width <- hb_obj_scurve2$a1
start_point <- bin_width/2

# Create bins and calculate the mean value for each bin range
bins <- seq(start_point, max(distance_scurve$distance) + bin_width,
            by = bin_width)
bin_labels <- bins[-length(bins)] + bin_width / 2  # mean of each bin range

distance_scurve <- distance_scurve |>
  mutate(bin = cut(distance, breaks = bins, include.lowest = TRUE, labels = bin_labels))

distance_hist <- ggplot(
  distance_scurve, aes(x = distance)) +
  geom_histogram(
    aes(y = after_stat(count) / sum(after_stat(count))),
    binwidth = bin_width) +
  scale_x_continuous(breaks = bin_labels,
                     labels = c(
                       expression(a[1]),
                       expression(2*a[1]),
                       expression(3*a[1]),
                       expression(4*a[1]),
                       expression(5*a[1]),
                       expression(6*a[1]),
                       expression(7*a[1]),
                       expression(8*a[1]),
                       expression(9*a[1]),
                       expression(10*a[1]),
                       expression(11*a[1]),
                       expression(12*a[1]),
                       expression(13*a[1])
                     )) +
  # geom_text(
  #   data = tibble(x = 12*bin_width, y = 0.7,
  #                 text = expression(a[1])),
  #   aes(
  #     x=x,
  #     y=y,
  #     label = text
  #   )
  #
  # ) +
  geom_vline(xintercept = start_point + 2*bin_width,
             linetype="solid",
             color = "#bdbdbd",
             linewidth=1) +
  ylab("Proportion of the number of edges") +
  xlab(expression(d^{(2)})) +
  ggtitle("(d)") +
  theme_minimal() +
  theme(aspect.ratio = 0.75,
        panel.border = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 12, hjust = 0.5, vjust = -0.5),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_line())

# distance_scurve_summary <- distance_scurve |>
#   group_by(bin) |>
#   summarise(count = n()) |>
#   mutate(std_counts = count/sum(count))

```

```{r}
#| fig-width: 6
#| fig-height: 4
#| fig-cap: "long edge removal"

trimesh_scurve + trimesh_scurve_removed2 + trimesh_scurve_removed1 +
  plot_layout(ncol = 3)
```

