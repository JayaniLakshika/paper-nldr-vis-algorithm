library(dplyr)
library(langevitour)
library(crosstalk)
library(quollr)
library(plotly)
library(readr)
library(DT)

## Import data
training_data_mnist <- read_rds("data/mnist/mnist_10_pcs_of_digit_1.rds")
training_data_mnist <- training_data_mnist |>
  mutate(ID = 1:NROW(training_data_mnist))

data_mnist <- training_data_mnist |>
  select(-ID) |>
  mutate(type = "data")

tsne_mnist <- read_rds("data/mnist/mnist_tsne30.rds")

mnist_scaled_obj <- gen_scaled_data(
  data = tsne_mnist)
tsne_minst_scaled <- mnist_scaled_obj$scaled_nldr

num_bins_x_mnist <- 19
lim1 <- mnist_scaled_obj$lim1
lim2 <- mnist_scaled_obj$lim2
r2_mnist <- diff(lim2)/diff(lim1)

mnist_model <- fit_highd_model(
  training_data = training_data_mnist,
  emb_df = tsne_minst_scaled,
  bin1 = num_bins_x_mnist,
  r2 = r2_mnist,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "PC"
)

df_bin_centroids_mnist <- mnist_model$df_bin_centroids
df_bin_mnist <- mnist_model$df_bin

## Triangulate bin centroids
tr1_object_mnist <- tri_bin_centroids(
  df_bin_centroids_mnist, x = "c_x", y = "c_y")
tr_from_to_df_mnist <- gen_edges(
  tri_object = tr1_object_mnist)

## Compute 2D distances
distance_mnist <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_mnist,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_mnist <- find_lg_benchmark(
  distance_edges = distance_mnist,
  distance_col = "distance")

sc_ltr_pos_mnist <- c(0.96, 0.96)

tr_df <- distinct(tibble::tibble(
  x = c(tr_from_to_df_mnist[["x_from"]], tr_from_to_df_mnist[["x_to"]]),
  y = c(tr_from_to_df_mnist[["y_from"]], tr_from_to_df_mnist[["y_to"]])))

distance_df_small_edges_mnist <- distance_mnist |>
  filter(distance < benchmark_mnist)

tr_from_to_df_mnist <- inner_join(
  tr_from_to_df_mnist, distance_df_small_edges_mnist,
  by = c("from", "to"))

trimesh_removed_mnist <- ggplot() +
  geom_segment(data = tr_from_to_df_mnist,
               aes(
                 x = x_from,
                 y = y_from,
                 xend = x_to,
                 yend = y_to),
               colour = "#33a02c",
               linewidth = 1) +
  geom_point(data = tsne_minst_scaled,
             aes(
               x = tSNE1,
               y = tSNE2
             ),
             alpha=0.1) +
  #interior_annotation("a", sc_ltr_pos_mnist) +
  theme(aspect.ratio = 1)

## hexagon binning to have regular hexagons
hb_obj_mnist <- hex_binning(
  data = tsne_minst_scaled,
  bin1 = num_bins_x_mnist,
  r2 = r2_mnist,
  q = 0.1)

tsne_data_with_hb_id <- hb_obj_mnist$data_hb_id

df_all_mnist <- dplyr::bind_cols(training_data_mnist |> select(-ID),
                                  tsne_data_with_hb_id)

df_b_mnist <- df_bin_mnist |>
  dplyr::filter(hb_id %in% df_bin_centroids_mnist$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b_mnist <- df_b_mnist[match(df_bin_centroids_mnist$hexID, df_b_mnist$hb_id),] |>
  dplyr::select(-hb_id)

# Apply the scaling
df_model_data_mnist <- bind_rows(data_mnist, df_b_mnist)
scaled_mnist <- scale_data_manual(df_model_data_mnist, "type") |>
  as_tibble()

scaled_mnist_data <- scaled_mnist |>
  filter(type == "data") |>
  select(-type)

scaled_mnist_data <- bind_cols(scaled_mnist_data,
                               tsne_minst_scaled |> select(-ID))

scaled_mnist_data_model <- scaled_mnist |>
  filter(type == "model") |>
  select(-type)

model_2d <- df_bin_centroids_mnist |>
  select(c_x, c_y) |>
  rename(c("tSNE1" = "c_x",
           "tSNE2" = "c_y"))

scaled_mnist_data_model <- bind_cols(scaled_mnist_data_model, model_2d)


# Combine with the true model for visualization
df <- dplyr::bind_rows(scaled_mnist_data_model |> mutate(type = "model"),
                       scaled_mnist_data |> mutate(type = "data"))


# df_all_mnist_n <- df_all_mnist |>
#   select(-ID, -hb_id) |>
#   dplyr::mutate(type = "data")

# df_b <- df_bin_scurve |>
#   dplyr::filter(hb_id %in% df_bin_centroids_mnist$hexID)
#
# ## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
# df_b <- df_b[match(df_bin_centroids_scurve$hexID, df_b$hb_id),] |>
#   dplyr::select(-hb_id) |>
#   dplyr::mutate(type = "model")
#
# df_exe <- dplyr::bind_rows(df_b, df_all_scurve_n)

## Set the maximum difference as the criteria
distance_df_small_edges_mnist <- distance_mnist |>
  dplyr::filter(distance < benchmark_mnist)

shared_df_mnist <- SharedData$new(df)

nldr_scurve <- shared_df_mnist |>
  ggplot(aes(x = tSNE1, y = tSNE2)) +
  geom_point(alpha=0.5, colour="#000000", size = 0.5) +
  theme_linedraw() +
  theme(
    #aspect.ratio = 1,
    plot.background = element_rect(fill = 'transparent', colour = NA),
    plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
    panel.background = element_rect(fill = 'transparent',
                                    colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )

nldr_scurve_plt <- ggplotly(nldr_scurve, width = as.character(round(600/r2_mnist, 0)),
                            height = "600", tooltip = "none") |>
  style(unselected=list(marker=list(opacity=1))) |>
  highlight(on="plotly_selected", off="plotly_deselect") #|>
  #config(displayModeBar = FALSE)


# langevitour_output <- langevitour::langevitour(df_exe[1:7],
#                                                lineFrom = distance_df_small_edges$from,
#                                                lineTo = distance_df_small_edges$to,
#                                                group = df_exe$type, pointSize = append(rep(2, NROW(df_b)), rep(1, NROW(df))),
#                                                levelColors = c("#000000", "#33a02c"),
#                                                link=shared_df_scurve,
#                                                link_filter=FALSE)

langevitour_output <- langevitour(df |> dplyr::select(starts_with("PC")),
                                  lineFrom = distance_df_small_edges_mnist$from,
                                  lineTo = distance_df_small_edges_mnist$to,
                                  group = df$type,
                                  pointSize = append(rep(1.5, NROW(scaled_mnist_data_model)), rep(1, NROW(scaled_mnist_data))),
                                  levelColors = c("#000000", "#33a02c"),
                                  lineColors = rep("#33a02c", nrow(distance_df_small_edges_mnist)),
                                  link=shared_df_mnist,
                                  link_filter=FALSE)

# Create a table widget
# datatableWidget <- datatable(
#   shared_df_mnist,
#   rownames=FALSE, width="100%",
#   class='compact cell-border hover', extensions='Buttons',
#   options=list(dom='Bfrtip',buttons=c('copy','csv','excel')))


linked_plt <- bscols(
  nldr_scurve_plt,
  langevitour_output,
  #datatableWidget,
  widths = c(6, 6),
  device = "sm"
)

linked_plt
