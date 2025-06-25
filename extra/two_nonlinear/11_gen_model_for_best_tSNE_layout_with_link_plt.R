## Import necessary libraries
library(quollr)
library(dplyr)
library(readr)
library(langevitour)
library(conflicted)
library(ggplot2)
library(colorspace)

library(plotly)
library(crosstalk)

conflicts_prefer(dplyr::filter)

## Import data
training_data_two_nonlinear_clusters <- read_rds("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_data.rds")

# training_data_two_nonlinear_clusters <- training_data_two_nonlinear_clusters |>
#   mutate(across(everything(), ~ (. - mean(.)) / sd(.)))

training_data_two_nonlinear_clusters <- training_data_two_nonlinear_clusters |>
  mutate(ID = 1:NROW(training_data_two_nonlinear_clusters))

tSNE_two_nonlinear_clusters <- read_rds("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_30.rds")

two_nonlinear_clusters_scaled_obj <- gen_scaled_data(
  data = tSNE_two_nonlinear_clusters)
tSNE_two_nonlinear_clusters_scaled <- two_nonlinear_clusters_scaled_obj$scaled_nldr |>
  mutate(ID = 1:NROW(tSNE_two_nonlinear_clusters))

## Compute hexbin parameters
num_bins_x_two_nonlinear_clusters <- 27
lim1 <- two_nonlinear_clusters_scaled_obj$lim1
lim2 <- two_nonlinear_clusters_scaled_obj$lim2
r2_two_nonlinear_clusters <- diff(lim2)/diff(lim1)

## Hexagonal binning to have regular hexagons
hb_obj_two_nonlinear_clusters <- hex_binning(
  data = tSNE_two_nonlinear_clusters_scaled,
  bin1 = num_bins_x_two_nonlinear_clusters,
  r2 = r2_two_nonlinear_clusters)

bin_width <- hb_obj_two_nonlinear_clusters$a1

all_centroids_df <- hb_obj_two_nonlinear_clusters$centroids
counts_df <- hb_obj_two_nonlinear_clusters$std_cts
tSNE_data_with_hb_id <- hb_obj_two_nonlinear_clusters$data_hb_id

df_bin_centroids_two_nonlinear_clusters <- extract_hexbin_centroids(centroids_df = all_centroids_df,
                                                  counts_df = counts_df) |>
  filter(drop_empty == FALSE)

df_all_two_nonlinear_clusters <- dplyr::bind_cols(training_data_two_nonlinear_clusters |> dplyr::select(-ID),
                                                  tSNE_data_with_hb_id)

df_bin_two_nonlinear_clusters <- avg_highd_data(data = df_all_two_nonlinear_clusters, col_start = "x")

## Triangulate bin centroids
tr1_object_two_nonlinear_clusters <- tri_bin_centroids(
  df_bin_centroids_two_nonlinear_clusters, x = "c_x", y = "c_y")
tr_from_to_df_two_nonlinear_clusters <- gen_edges(
  tri_object = tr1_object_two_nonlinear_clusters)

## Compute 2D distances
distance_two_nonlinear_clusters <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_two_nonlinear_clusters,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_two_nonlinear_clusters <- find_lg_benchmark(
  distance_edges = distance_two_nonlinear_clusters,
  distance_col = "distance")

tr_df <- distinct(tibble::tibble(
  x = c(tr_from_to_df_two_nonlinear_clusters[["x_from"]], tr_from_to_df_two_nonlinear_clusters[["x_to"]]),
  y = c(tr_from_to_df_two_nonlinear_clusters[["y_from"]], tr_from_to_df_two_nonlinear_clusters[["y_to"]])))

distance_df_small_edges_two_nonlinear_clusters <- distance_two_nonlinear_clusters |>
  filter(distance < benchmark_two_nonlinear_clusters)

tr_from_to_df_two_nonlinear_clusters <- inner_join(
  tr_from_to_df_two_nonlinear_clusters, distance_df_small_edges_two_nonlinear_clusters,
  by = c("from", "to"))

### Define type column
df <- df_all_two_nonlinear_clusters |>
  dplyr::select(tidyselect::starts_with("x")) |>
  dplyr::mutate(type = "data") ## original dataset

df_b <- df_bin_two_nonlinear_clusters |>
  dplyr::filter(hb_id %in% df_bin_centroids_two_nonlinear_clusters$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_two_nonlinear_clusters$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

df_exe <- dplyr::bind_rows(df_b, df)

## Set the maximum difference as the criteria
# distance_df_small_edges_two_nonlinear_clusters <- distance_two_nonlinear_clusters |>
#   dplyr::filter(distance < benchmark_two_nonlinear_clusters)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges_two_nonlinear_clusters$from,
                         lineTo = distance_df_small_edges_two_nonlinear_clusters$to,
                         group = df_exe$type, pointSize = append(rep(1, NROW(df_b)), rep(0.5, NROW(df))),
                         levelColors = c("#6a3d9a", "#33a02c"))

#### With scaled data

data_two_nonlinear_clusters <- training_data_two_nonlinear_clusters |>
  select(-ID) |>
  mutate(type = "data")

# Apply the scaling
df_model_data <- bind_rows(data_two_nonlinear_clusters, df_b)
scaled_two_nonlinear_clusters <- scale_data_manual(df_model_data, "type") |>
  as_tibble()

scaled_two_nonlinear_clusters_data <- scaled_two_nonlinear_clusters |>
  filter(type == "data") |>
  select(-type)

scaled_two_nonlinear_clusters_data_model <- scaled_two_nonlinear_clusters |>
  filter(type == "model") |>
  select(-type)


df_b_two_nonlinear_clusters <- df_bin_two_nonlinear_clusters |>
  dplyr::filter(hb_id %in% df_bin_centroids_two_nonlinear_clusters$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b_two_nonlinear_clusters <- df_b_two_nonlinear_clusters[match(df_bin_centroids_two_nonlinear_clusters$hexID, df_b_two_nonlinear_clusters$hb_id),] |>
  dplyr::select(-hb_id) |>
  select(-type)

# Combine with the true model for visualization
df <- scaled_two_nonlinear_clusters_data |> mutate(type = "data") |>
  bind_cols(tSNE_data_with_hb_id) |>
  select(-ID, -hb_id) |>
  dplyr::bind_rows(scaled_two_nonlinear_clusters_data_model |> mutate(type = "model"))

shared_df <- crosstalk::SharedData$new(df)

## Set the maximum difference as the criteria
distance_df_small_edges_two_nonlinear_clusters <- distance_two_nonlinear_clusters |>
  dplyr::filter(distance < benchmark_two_nonlinear_clusters)

# Visualize with langevitour
langevitour_output <- langevitour(df |> dplyr::select(tidyselect::starts_with("x")),
            lineFrom = distance_df_small_edges_two_nonlinear_clusters$from,
            lineTo = distance_df_small_edges_two_nonlinear_clusters$to,
            group = df$type,
            pointSize = append(rep(1.5, NROW(scaled_two_nonlinear_clusters_data_model)), rep(1, NROW(scaled_two_nonlinear_clusters_data))),
            levelColors = c("#000000", "#33a02c"),
            lineColors = rep("#33a02c", nrow(distance_df_small_edges_two_nonlinear_clusters)),
            link=shared_df,
            link_filter=FALSE)


nldr_plt <- shared_df |>
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

nldr_plt <- ggplotly(nldr_plt, width = as.character(round(600/r2_two_nonlinear_clusters, 0)),
                     height = "600", tooltip = "none") |>
  style(unselected=list(marker=list(opacity=1))) |>
  highlight(on="plotly_selected", off="plotly_deselect") |>
  config(displayModeBar = FALSE)

linked_plt <- crosstalk::bscols(
  nldr_plt,
  langevitour_output,
  widths = c(6, 6),
  device = "sm"
)

linked_plt
