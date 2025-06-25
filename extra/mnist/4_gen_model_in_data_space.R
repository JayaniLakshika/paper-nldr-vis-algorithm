## Import necessary libraries
library(quollr)
library(dplyr)
library(readr)
library(langevitour)

## Import data
training_data_mnist <- read_rds("data/mnist/mnist_10_pcs_of_digit_1.rds")
training_data_mnist <- training_data_mnist |>
  mutate(ID = 1:NROW(training_data_mnist))

tsne_minst <- read_rds("data/mnist/mnist_tsne89.rds")

mnist_scaled_obj <- gen_scaled_data(
  data = tsne_minst)
tsne_minst_scaled <- mnist_scaled_obj$scaled_nldr

## Compute hexbin parameters
num_bins_x_mnist <- 30
lim1 <- mnist_scaled_obj$lim1
lim2 <- mnist_scaled_obj$lim2
r2_mnist <- diff(lim2)/diff(lim1)

## Hexagonal binning to have regular hexagons
hb_obj_mnist <- hex_binning(
  data = tsne_minst_scaled,
  bin1 = num_bins_x_mnist,
  r2 = r2_mnist)

all_centroids_df <- hb_obj_mnist$centroids
counts_df <- hb_obj_mnist$std_cts
tsne_data_with_hb_id <- hb_obj_mnist$data_hb_id

df_bin_centroids_mnist <- extract_hexbin_centroids(centroids_df = all_centroids_df,
                                             counts_df = counts_df) |>
  filter(drop_empty == FALSE)

## first quartile used as the default
benchmark_to_rm_lwd_hex <- quantile(df_bin_centroids_mnist$std_counts,
                                    probs = c(0,0.25,0.5,0.75,1))[2]

## To identify low density hexagons
df_bin_centroids_low <- df_bin_centroids_mnist |>
  filter(std_counts <= benchmark_to_rm_lwd_hex)

## To identify low-density hexagons needed to remove by investigating neighbouring mean density
identify_rm_bins <- find_low_dens_hex(df_bin_centroids_all = df_bin_centroids_mnist,
                                      bin1 = num_bins_x_mnist,
                                      df_bin_centroids_low = df_bin_centroids_low)

## To remove low-density hexagons
df_bin_centroids_mnist <- df_bin_centroids_mnist |>
  filter(!(hexID %in% identify_rm_bins))

df_all_mnist <- dplyr::bind_cols(training_data_mnist |> dplyr::select(-ID),
                                 tsne_data_with_hb_id)

df_bin_mnist <- avg_highd_data(data = df_all_mnist, col_start = "PC")

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


### Define type column
df <- df_all_mnist |>
  dplyr::select(tidyselect::starts_with("PC")) |>
  dplyr::mutate(type = "data") ## original dataset

df_b <- df_bin_mnist |>
  dplyr::filter(hb_id %in% df_bin_centroids_mnist$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_mnist$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

df_exe <- dplyr::bind_rows(df_b, df)

## Set the maximum difference as the criteria
distance_df_small_edges_mnist <- distance_mnist |>
  dplyr::filter(distance < benchmark_mnist)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges_mnist$from,
                         lineTo = distance_df_small_edges_mnist$to,
                         group = df_exe$type, pointSize = append(rep(1, NROW(df_b)), rep(0.5, NROW(df))),
                         levelColors = c("#6a3d9a", "#33a02c"))

#### With scaled data

data_mnist <- training_data_mnist |>
  bind_cols(tsne_minst_scaled |> select(-ID)) |>
  select(-ID)

data_mnist <- data_mnist |>
  mutate(type = if_else((tSNE1 <= 1) & (tSNE1 >= 0.78) & (tSNE2 <= 0.62) & (tSNE2 >= 0.39), "small_clust", "big_clust"))

data_mnist <- data_mnist |>
  select(-tSNE1, -tSNE2)

# Apply the scaling
df_model_data <- bind_rows(data_mnist, df_b)
scaled_mnist <- scale_data_manual(df_model_data, "type") |>
  as_tibble()

scaled_mnist_data <- scaled_mnist |>
  filter(type %in% c("small_clust", "big_clust")) #|>
  #select(-type)

scaled_mnist_data_model <- scaled_mnist |>
  filter(type == "model") #|>
  #select(-type)


df_b_mnist <- df_bin_mnist |>
  dplyr::filter(hb_id %in% df_bin_centroids_mnist$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b_mnist <- df_b_mnist[match(df_bin_centroids_mnist$hexID, df_b_mnist$hb_id),] |>
  dplyr::select(-hb_id) |>
  select(-type)

# Combine with the true model for visualization
df <- dplyr::bind_rows(scaled_mnist_data_model,
                       scaled_mnist_data)

## Set the maximum difference as the criteria
distance_df_small_edges_mnist <- distance_mnist |>
  dplyr::filter(distance < benchmark_mnist)

# Visualize with langevitour
langevitour(df |> dplyr::select(-type),
            lineFrom = distance_df_small_edges_mnist$from,
            lineTo = distance_df_small_edges_mnist$to,
            group = df$type,
            pointSize = append(rep(1.5, NROW(scaled_mnist_data_model)), rep(1, NROW(scaled_mnist_data))),
            levelColors = c("#999999","#000000", '#ff7f00'),
            lineColors = rep("#000000", nrow(distance_df_small_edges_mnist)))


