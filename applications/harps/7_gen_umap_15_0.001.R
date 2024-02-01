library(readr)
library(dplyr)
library(umap)
library(ggplot2)

set.seed(20240110)
source("quollr_code.R", local = TRUE)

## Import data
training_data <- read_rds("data/harps/harps_data.rds")
training_data <- training_data |>
  mutate(ID = 1:NROW(training_data))

clusters <- training_data$cluster
select_names <-  names(training_data)[-c(18, length(training_data))]

training_data <- training_data |> select(-cluster)
names(training_data) <- append(paste0("x", 1:17), "ID")

## Run only once

UMAP_fit <- umap(training_data |> dplyr::select(-ID), n_neighbors = 15, n_components =  2, min_dist = 0.001)

UMAP_data <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data)[1:(ncol(UMAP_data))] <- paste0(rep("UMAP",(ncol(UMAP_data))), 1:(ncol(UMAP_data)))

UMAP_data_harps <- UMAP_data |>
  mutate(ID = training_data$ID)

write_rds(UMAP_data_harps, "data/harps/harps_umap_15_0.001.rds")

UMAP_data_harps <- read_rds("data/harps/harps_umap_15_0.001.rds")

UMAP_data_harps |>
  ggplot(aes(x = UMAP2,
             y = UMAP1,
             color = clusters)) +
  geom_point() #+
#geom_point(aes(shape=clusters, color=clusters)) +
#scale_shape_manual(values=c(16, 5, 5, 13, 8, 15, 16, 15, 10, 17, 24, 25, 17, 17)) +  ggtitle("Visualization from generated \n results by UMAP") + xlab("t-SNE X dimension") + ylab("t-SNE Y dimension") +
#theme(legend.position = "none") + coord_fixed()

# UMAP_data_harps <- UMAP_data_harps |>
#   select(UMAP1, UMAP2) |>
#   mutate(ID = 1:NROW(UMAP_data_harps))

### UMAP

num_bins_umap_harps <- calculate_effective_x_bins(.data = UMAP_data_harps, x = UMAP1,
                                                  cell_area = 1) #58

num_bins_umap_harps <- 43


shape_val_umap_harps <- calculate_effective_shape_value(.data = UMAP_data_harps,
                                                        x = UMAP1, y = UMAP2)

shape_val_tsne_harps <- 0.4466309

## To extract bin centroids
hexbin_data_object_umap_harps <- extract_hexbin_centroids(nldr_df = UMAP_data_harps, num_bins = num_bins_umap_harps,
                                                          shape_val = shape_val_umap_harps, x = UMAP1, y = UMAP2)

df_bin_centroids_umap_harps <- hexbin_data_object_umap_harps$hexdf_data

## Identify bins with low-density
# identify_rm_bins <- find_low_density_hexagons(df_bin_centroids_umap_harps, num_bins_umap_harps, benchmark_rm_hex = 0.01)
#
# df_bin_centroids_umap_harps <- df_bin_centroids_umap_harps |>
#   filter(!(hexID %in% identify_rm_bins))

UMAP_data_with_hb_id_harps <- UMAP_data_harps |>
  dplyr::mutate(hb_id = hexbin_data_object_umap_harps$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_umap_harps <- dplyr::bind_cols(training_data |> dplyr::select(-ID), UMAP_data_with_hb_id_harps)

## Averaged on high-D
df_bin_umap_harps <- avg_highD_data(.data = df_all_umap_harps)

## Triangulate bin centroids
tr1_object_umap_harps <- triangulate_bin_centroids(df_bin_centroids_umap_harps, x, y)
tr_from_to_df_umap_harps <- generate_edge_info(triangular_object = tr1_object_umap_harps)

## Compute 2D distances
distance_umap_harps <- cal_2D_dist(.data = tr_from_to_df_umap_harps)

## To find the benchmark value
benchmark_umap_harps <- find_benchmark_value(.data = distance_umap_harps, distance_col = distance)
#benchmark_umap_harps <- 1.900017


trimesh_removed_umap_harps <- remove_long_edges(.data = distance_umap_harps, benchmark_value = benchmark_umap_harps,
                                                triangular_object = tr1_object_umap_harps, distance_col = distance)

trimesh_removed_umap_harps <- trimesh_removed_umap_harps +
  # xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  # theme(axis.text = element_text(size = 5),
  #       axis.title = element_text(size = 7)) +
  geom_point(colour = "#33a02c", size = 0.5) +
  theme_linedraw() +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate(geom = 'text', label = 'a', x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 3)

#select_names <- training_data_1 |> dplyr::select(-c(ID, cluster)) |> names()

############ Compute absolute residuals

df_bin_train <- df_bin_umap_harps
names(df_bin_train)[-1] <- paste0("avg_", names(df_bin_train)[-1])

error_df <- df_all_umap_harps |>
  dplyr::left_join(df_bin_train, by = c("hb_id" = "hb_id")) ## Map high-D averaged/weighted mean coordinates

# prediction_df_join <- prediction_df_join |>
#   dplyr::left_join(data, by = c("ID" = "ID")) ## Map high-D data

for (i in 1:(NCOL(df_bin_train) - 1)) {

  error_df[ , paste0("abs_residual_", "x", i)] <- abs(error_df[ , paste0("x", i)] - error_df[ , paste0("avg_", "x", i)])

}

error_df <- error_df |>
  dplyr::mutate(total = rowSums(dplyr::pick(tidyselect::starts_with(paste0("abs_residual_", "x")))))

library(ggbeeswarm)
error_df$group <- "1"
ggplot(error_df, aes(x = group, y = total)) +
  geom_quasirandom()+
  ylim(0, max(unlist(error_df$total))+ 0.5) + coord_flip()

### The minimum error is 0 and the maximum is 1.38
### There is lot of points with error 0,

error_df <- error_df |>
  mutate(type = if_else(total <= 0, "no error",
                        if_else(total <= 0.5, "error 0-0.5",
                                if_else(total <= 1, "error 0.5-1",
                                        if_else(total <= 1.5, "error 1-1.5",
                                                if_else(total <= 2, "error 1.5-2",
                                                        "error greter than 2"))))))

error_sum_df <- tibble::tibble(n_neighbors = 15, min_dist = 0.001, metric = "euclidean", total_error = sum(error_df$total))
write_csv(error_sum_df, "data/harps/error_harps_umap.csv", append = TRUE)



####
### Define type column
df <- df_all_umap_harps |>
  dplyr::select(tidyselect::starts_with("x")) #|>
#dplyr::mutate(type = "data")  ## original dataset

residual_df <- error_df |> select(type)

df <- dplyr::bind_cols(df, residual_df)

df_b <- df_bin_umap_harps |>
  dplyr::filter(hb_id %in% df_bin_centroids_umap_harps$hexID) |>
  dplyr::select(-hb_id) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

df_exe <- dplyr::bind_rows(df_b, df)
names(df_exe) <- append(select_names, "type")

distance_df_small_edges <- distance_umap_harps %>%
  dplyr::filter(distance < benchmark_umap_harps)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to, group = factor(df_exe$type, levels = c("no error", "error 0-0.5", "error 0.5-1", "error 1-1.5", "error 1.5-2", "error greter than 2", "model")), pointSize = 3,
                         levelColors = c("#b15928", "#1f78b4",
                                         "#fb9a99", "#cab2d6", "#ff7f00",
                                         "#e31a1c", "#33a02c"))



#tour_umap_harps <- show_langevitour(df_all_umap_harps, df_bin_umap_harps, df_bin_centroids_umap_harps, benchmark_value = benchmark_umap_harps, distance = distance_umap_harps, distance_col = distance)
