library(readr)

### Spin-brush
#remotes::install_github("casperhart/detourr")
library(detourr)

training_data_pbmc <- read_rds("data/pbmc/pbmc_3k_festem/pbmc_pca_50.rds")
training_data_pbmc <- training_data_pbmc[, 1:15] |>
  mutate(ID = 1:NROW(training_data_pbmc))

detour(training_data_pbmc,
       tour_aes(projection = PC_1:PC_15)) |>
  tour_path(grand_tour(2), fps = 60,
            max_bases=20) |>
  show_scatter(alpha = 0.7,
               axes = FALSE)

## Import detour results
detour_results1 <- read_csv("data/pbmc/pbmc_3k_festem/detourr_export.csv")
#names(detour_results)[1:15] <- paste0("PC_", 1:15)
table(training_data_pbmc$cell_label, detour_results1$colour)

## Import detour results
detour_results2 <- read_csv("data/pbmc/pbmc_3k_festem/detourr_export2.csv")
table(training_data_pbmc$cell_label, detour_results2$colour)

## Import detour results
detour_results3 <- read_csv("data/pbmc/pbmc_3k_festem/detourr_export3.csv")
table(training_data_pbmc$cell_label, detour_results3$colour)


### Spin-brush with detour and model
UMAP_pbmc <- read_rds("data/pbmc/pbmc_umap_5_min_dist_0.99_metric_cosine.rds")
shape_value_pbmc <- calculate_effective_shape_value(.data = UMAP_pbmc,
                                                    x = UMAP1, y = UMAP2) ## 0.8772751

num_non_empty_bins_pbmc <- 83
num_bins_pbmc <- find_non_empty_bins(nldr_df = UMAP_pbmc, x = "UMAP1", y = "UMAP2", shape_val = shape_value_pbmc,
                                     non_empty_bins = num_non_empty_bins_pbmc)

## To extract bin centroids
hexbin_data_object_pbmc <- extract_hexbin_centroids(nldr_df = UMAP_pbmc,
                                                    num_bins = num_bins_pbmc,
                                                    shape_val = shape_value_pbmc, x = UMAP1, y = UMAP2)

df_bin_centroids_pbmc <- hexbin_data_object_pbmc$hexdf_data

UMAP_pbmc_with_hb_id <- UMAP_pbmc |>
  dplyr::mutate(hb_id = hexbin_data_object_pbmc$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_pbmc <- dplyr::bind_cols(training_data_pbmc |> dplyr::select(-ID), UMAP_pbmc_with_hb_id)

## Averaged on high-D
df_bin_pbmc <- avg_highD_data(.data = df_all_pbmc, column_start_text = "PC")

## Triangulate bin centroids
tr1_object_pbmc <- triangulate_bin_centroids(df_bin_centroids_pbmc, x, y)
tr_from_to_df_pbmc <- generate_edge_info(triangular_object = tr1_object_pbmc)

## Compute 2D distances
distance_pbmc <- cal_2d_dist(.data = tr_from_to_df_pbmc)

## To find the benchmark value
benchmark_pbmc <- find_benchmark_value(.data = distance_pbmc, distance_col = distance)

df_bin_train <- df_bin_pbmc
names(df_bin_train)[-1] <- paste0("avg_", names(df_bin_train)[-1])

error_df <- df_all_pbmc |>
  dplyr::left_join(df_bin_train, by = c("hb_id" = "hb_id")) ## Map high-D averaged/weighted mean coordinates

# prediction_df_join <- prediction_df_join |>
#   dplyr::left_join(data, by = c("ID" = "ID")) ## Map high-D data

for (i in 1:(NCOL(df_bin_train) - 1)) {

  error_df[ , paste0("abs_residual_", "PC_", i)] <- abs(error_df[ , paste0("PC_", i)] - error_df[ , paste0("avg_", "PC_", i)])

}

error_df <- error_df |>
  dplyr::mutate(total = rowSums(dplyr::pick(tidyselect::starts_with(paste0("abs_residual_", "PC_")))))

# library(ggbeeswarm)
# error_df$group <- "1"
# ggplot(error_df, aes(x = group, y = total)) +
#   geom_quasirandom()+
#   ylim(0, max(unlist(error_df$total))+ 0.5) + coord_flip()

### The minimum error is 0 and the maximum is 42.17439
### There is lot of points with error 0,

error_df <- error_df |>
  mutate(type = if_else(total <= 2, "error 2 or less",
                        if_else(total <= 10, "error 2-10",
                                if_else(total <= 15, "error 10-15",
                                        if_else(total <= 20, "error 15-20",
                                                if_else(total <= 25, "error 20-25",
                                                        if_else(total <= 30, "error 25-30",
                                                                if_else(total <= 35, "error 30-35", "error greter than 35"))))))))

### Define type column
df <- df_all_pbmc |>
  dplyr::select(tidyselect::starts_with("PC")) #|>
#dplyr::rename("type" = "cell_label") ## original dataset

residual_df <- error_df |> select(type)

df <- dplyr::bind_cols(df, residual_df)

df_b <- df_bin_pbmc |>
  dplyr::filter(hb_id %in% df_bin_centroids_pbmc$hexID) |>
  dplyr::select(-hb_id) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

df_exe <- dplyr::bind_rows(df_b, df)

distance_df_small_edges <- distance_pbmc %>%
  dplyr::filter(distance < benchmark_pbmc)

detour(df_exe,
       tour_aes(projection = PC_1:PC_15,
                colour = type)) |>
  tour_path(grand_tour(2), fps = 60,
            max_bases=20) |>
  show_scatter(alpha = 0.7,
               axes = FALSE,
               edges = distance_df_small_edges |> select(from, to) |> as.matrix(),
               palette = c("#fee0d2", "#fcbba1",
                           "#fc9272", "#fb6a4a", "#ef3b2c",
                           "#cb181d", "#a50f15", "#99000d", "#33a02c"))


### Experiment with liminal
library(liminal)

training_data_pbmc <- training_data_pbmc |>
  mutate(cell_label = UMAP_pbmc$cell_label)

limn_tour_link(embed_data = UMAP_pbmc |> select("UMAP1", "UMAP2"),
               tour_data = training_data_pbmc,
               cols = PC_1:PC_15, # tour columns to select
               color = cell_label # variable to highlight across both view, can come for either data frames
)
