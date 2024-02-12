library(readr)

set.seed(20240110)

source("quollr_code.R", local = TRUE)
source("nldr_code.R", local = TRUE)

## Import data
training_data <- read_rds("data/sphere/sphere_training.rds")
UMAP_data_sphere <- read_rds("data/sphere/sphere_umap.rds")

### tSNE
tSNE_data_sphere <- Fit_tSNE(sphere_df |> distinct(), opt_perplexity = calculate_effective_perplexity(sphere_df |> distinct()),
                          with_seed = 20240110)

tSNE_data_sphere <- tSNE_data_sphere |>
  select(-ID) |>
  mutate(ID = training_data_5$ID)

plot_tSNE_2D(tSNE_data_sphere)

### UMAP

UMAP_fit <- umap(sphere_df |> distinct(), n_neighbors = 100, n_components =  2)

UMAP_data_sphere <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data_sphere)[1:(ncol(UMAP_data_sphere))] <- paste0(rep("UMAP",(ncol(UMAP_data_sphere))), 1:(ncol(UMAP_data_sphere)))

UMAP_data_sphere <- UMAP_data_sphere |>
  mutate(ID = training_data_5$ID)

plot_UMAP_2D(UMAP_data_sphere)

### Phate

PHATE_data <- Fit_PHATE(sphere_df |> distinct(), knn = 5, with_seed = 20240110)
PHATE_data <- PHATE_data |>
  select(PHATE1, PHATE2)
PHATE_data <- PHATE_data |>
  mutate(ID = training_data_5$ID)

plot_PHATE_2D(PHATE_data)

### TriMAP

tem_dir <- tempdir()

Fit_TriMAP_data(sphere_df |> distinct(), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_TriMAP_values.csv")

Fit_TriMAP(as.integer(2), as.integer(5), as.integer(4), as.integer(3), path, path2)

TriMAP_data <- read_csv(path2)
TriMAP_data <- TriMAP_data |>
  mutate(ID = training_data_5$ID)

plot_TriMAP_2D(TriMAP_data)

### PaCMAP

tem_dir <- tempdir()

Fit_PacMAP_data(sphere_df |> distinct(), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_PaCMAP_values.csv")

Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.9, as.integer(2), path, path2)

PacMAP_data <- read_csv(path2)
PacMAP_data <- PacMAP_data |>
  mutate(ID = training_data_5$ID)

plot_PaCMAP_2D(PacMAP_data)

## UMAP

num_bins_UMAP_data_sphere <- calculate_effective_x_bins(.data = UMAP_data_sphere, x = UMAP1, cell_area = 1)
shape_val_UMAP_data_sphere <- calculate_effective_shape_value(.data = UMAP_data_sphere,
                                                          x = UMAP1, y = UMAP2) ## 1.259938
## To extract bin centroids
hexbin_data_object_UMAP_data_sphere <- extract_hexbin_centroids(nldr_df = UMAP_data_sphere, num_bins = num_bins_UMAP_data_sphere, shape_val = shape_val_UMAP_data_sphere, x = UMAP1, y = UMAP2)

df_bin_centroids_UMAP_data_sphere <- hexbin_data_object_UMAP_data_sphere$hexdf_data

UMAP_data_with_hb_id_sphere <- UMAP_data_sphere |>
  dplyr::mutate(hb_id = hexbin_data_object_UMAP_data_sphere$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_UMAP_data_sphere <- dplyr::bind_cols(sphere_df |> distinct(), UMAP_data_with_hb_id_sphere)

## Averaged on high-D
df_bin_UMAP_data_sphere <- avg_highD_data(.data = df_all_UMAP_data_sphere)

## Triangulate bin centroids
tr1_object_UMAP_data_sphere <- triangulate_bin_centroids(df_bin_centroids_UMAP_data_sphere, x, y)
tr_from_to_df_UMAP_data_sphere <- generate_edge_info(triangular_object = tr1_object_UMAP_data_sphere)

# ggplot(df_bin_centroids_UMAP_data_sphere, aes(x = x, y = y)) +
#   geom_segment(data = tr_from_to_df_UMAP_data_sphere, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
#   geom_point(size = 2, colour = "#33a02c") +
#   coord_equal()


## Compute 2D distances
distance_UMAP_data_sphere <- cal_2d_dist(.data = tr_from_to_df_UMAP_data_sphere)


## To find the benchmark value
benchmark_UMAP_data_sphere <- find_benchmark_value(.data = distance_UMAP_data_sphere, distance_col = distance)



# colour_long_edges(.data = distance_UMAP_data_sphere, benchmark_value = benchmark_UMAP_data_sphere,
#                   triangular_object = tr1_object_UMAP_data_sphere, distance_col = distance)



trimesh_UMAP_data_sphere <- ggplot(df_bin_centroids_UMAP_data_sphere, aes(x = x, y = y)) +
  geom_segment(data = tr_from_to_df_UMAP_data_sphere, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
  geom_point(size = 2, colour = "#33a02c") +
  coord_equal()

trimesh_UMAP_data_sphere <- trimesh_UMAP_data_sphere +
  #ggtitle("(a)") +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme_light() +
  theme(legend.position = "none", plot.title = element_text(size = 5, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()#change legend key width
  ) +
  annotate(geom = 'text', label = 'a', x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 3)
# theme(axis.text = element_text(size = 5),
#       axis.title = element_text(size = 7))

trimesh_gr_UMAP_data_sphere <- colour_long_edges(.data = distance_UMAP_data_sphere, benchmark_value = benchmark_UMAP_data_sphere, triangular_object = tr1_object_UMAP_data_sphere, distance_col = distance)

trimesh_gr_UMAP_data_sphere <- trimesh_gr_UMAP_data_sphere +
  geom_point(size = 2, colour = "#33a02c") +
  #ggtitle("(b)") +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme_light() +
  #coord_equal() +
  theme(legend.position = "none", plot.title = element_text(size = 5, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()#change legend key width
  ) +
  annotate(geom = 'text', label = 'b', x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 3)


trimesh_removed_UMAP_data_sphere <- remove_long_edges(.data = distance_UMAP_data_sphere, benchmark_value = benchmark_UMAP_data_sphere, triangular_object = tr1_object_UMAP_data_sphere, distance_col = distance)

trimesh_removed_UMAP_data_sphere <- trimesh_removed_UMAP_data_sphere +
  geom_point(size = 2, colour = "#33a02c") +
  #ggtitle("(b)") +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme_light() +
  #coord_equal() +
  theme(legend.position = "none", plot.title = element_text(size = 5, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()#change legend key width
  ) +
  annotate(geom = 'text', label = 'b', x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 3)

show_langevitour(sphere_df |> distinct(), df_bin_UMAP_data_sphere, df_bin_centroids_UMAP_data_sphere, benchmark_value = benchmark_UMAP_data_sphere,
                             distance_df = distance_UMAP_data_sphere, distance_col = distance)

############ Compute absolute residuals

df_bin_train <- df_bin_UMAP_data_sphere
names(df_bin_train)[-1] <- paste0("avg_", names(df_bin_train)[-1])

error_df <- df_all_UMAP_data_sphere |>
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
  ylim(0, max(unlist(error_df$total))+ 0.05) + coord_flip()

### The minimum error is 0 and the maximum is 42.17439
### There is lot of points with error 0,

error_df <- error_df |>
  mutate(type = if_else(total == 0, "no error",
                        if_else(total <= 0.15, "error 0-0.15",
                                if_else(total <= 0.3, "error 0.15-0.3",
                                        if_else(total <= 0.45, "error 0.3-0.45",
                                                if_else(total <= 0.6, "error 0.45-0.6",
                                                        if_else(total <= 0.75, "error 0.6-0.75",
                                                                if_else(total <= 0.9, "error 0.75-0.9", "error greter than 0.9"))))))))

# error_sum_df <- tibble::tibble(n_neighbors = 354, min_dist = 0.99, metric = "euclidean", total_error = sum(error_df$total))
# write_csv(error_sum_df, "data/pbmc/pbmc_3k_festem/error_UMAP_data_sphere_umap.csv", append = TRUE)

### Define type column
df <- df_all_UMAP_data_sphere |>
  dplyr::select(tidyselect::starts_with("x")) #|>
#dplyr::rename("type" = "cell_label") ## original dataset

residual_df <- error_df |> select(type)

df <- dplyr::bind_cols(df, residual_df)


#df$type <- as.factor(df$type)

#levels(df$type) <- c("Memory \nCD4 T", "Naive CD4 T", "CD14+ Mono",  "B", "CD8 T", "FCGR3A+ \n Mono", "NK", "M-MDSC\n-like", "CD27-CD+ \n Memory T", "DC")

df_b <- df_bin_UMAP_data_sphere |>
  dplyr::filter(hb_id %in% df_bin_centroids_UMAP_data_sphere$hexID) |>
  dplyr::select(-hb_id) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

df_exe <- dplyr::bind_rows(df_b, df)

distance_df_small_edges <- distance_UMAP_data_sphere %>%
  dplyr::filter(distance < benchmark_UMAP_data_sphere)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to, group = factor(df_exe$type , levels = c("no error", "error 0-0.15", "error 0.15-0.3", "error 0.3-0.45", "error 0.45-0.6", "error 0.6-0.75", "error 0.75-0.9", "error greter than 0.9", "model")
                         ), pointSize = 3,
                         levelColors = c("#fee0d2", "#fcbba1",
                                         "#fc9272", "#fb6a4a", "#ef3b2c",
                                         "#cb181d", "#a50f15", "#99000d", "#33a02c"))
