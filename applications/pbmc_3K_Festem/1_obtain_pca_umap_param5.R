library(dplyr)
library(rsample)
library(readr)
library(umap)

set.seed(20240110)

source("quollr_code.R", local = TRUE)
source("nldr_code.R", local = TRUE)

#save(gene.list,label.list,umap.list,plots.list,file = "./results/pbmc3k_clustering_UMAP.RData")

data_pca <- read_rds(file = "data/pbmc/pbmc_3k_festem/pbmc_pca_50.rds")

data_pca <- data_pca[,1:15]

training_data_pbmc <- data_pca |>
  mutate(ID = row_number())

## UMAP

UMAP_fit <- umap(training_data_pbmc |> dplyr::select(-ID), n_neighbors = 5, n_components =  2, metric = "cosine", min_dist = 0.1)

UMAP_pbmc <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_pbmc)[1:(ncol(UMAP_pbmc))] <- paste0(rep("UMAP",(ncol(UMAP_pbmc))), 1:(ncol(UMAP_pbmc)))

UMAP_pbmc <- UMAP_pbmc |>
  mutate(ID = training_data_pbmc$ID)

UMAP_pbmc_old <- read_rds("data/pbmc/pbmc_3k_festem/pbmc_umap.rds")

UMAP_pbmc <- UMAP_pbmc |>
  mutate(cell_label = UMAP_pbmc_old$cell_label)

## Run only once
write_rds(UMAP_pbmc, file = "data/pbmc/pbmc_3k_festem/pbmc_umap_5_min_dist_0.1.rds")

plot_UMAP_2D(UMAP_pbmc) + #ggtitle("(b)") +
  theme_linedraw() +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate(geom = 'text', label = 'b', x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 3)

UMAP_pbmc <- read_rds("data/pbmc/pbmc_3k_festem/pbmc_umap_5_min_dist_0.1.rds")

num_bins_pbmc <- calculate_effective_x_bins(.data = UMAP_pbmc, x = UMAP1,
                                            cell_area = 1) ##23

shape_value_pbmc <- calculate_effective_shape_value(.data = UMAP_pbmc,
                                                    x = UMAP1, y = UMAP2) ## 0.8772751

## To extract bin centroids
hexbin_data_object_pbmc <- extract_hexbin_centroids(nldr_df = UMAP_pbmc,
                                                    num_bins = num_bins_pbmc,
                                                    shape_val = shape_value_pbmc, x = UMAP1, y = UMAP2)

df_bin_centroids_pbmc <- hexbin_data_object_pbmc$hexdf_data

##########

## Data set with all possible centroids in the hexagonal grid

full_centroid_df <- generate_full_grid_centroids(df_bin_centroids_pbmc)

## Generate all coordinates of hexagons
hex_grid <- full_hex_grid(full_centroid_df)

hex_full_count_df <- generate_full_grid_info(df_bin_centroids_pbmc)

#write_rds(hex_full_count_df, "data/pbmc/pbmc_3k_festem/pbmc_umap_hex_23.rds")

##########

ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_point(data = UMAP_pbmc, aes(x = UMAP1, y = UMAP2), alpha = 0.5) +
  scale_fill_viridis_c(direction = -1, na.value = "#ffffff", option = "C") +
  xlim(-5, 7) +
  ylim(-10, 10) +
  theme_void() +
  theme(legend.position="none", legend.direction="horizontal", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6)) +
  guides(fill = guide_colourbar(title = "Standardized count")) +
  annotate(geom = 'text', label = "a", x = -Inf, y = Inf, hjust = -0.3, vjust = 1, size = 3)


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
distance_pbmc <- cal_2D_dist(.data = tr_from_to_df_pbmc)

## To find the benchmark value
benchmark_pbmc <- find_benchmark_value(.data = distance_pbmc, distance_col = distance)
benchmark_pbmc <- 1.055271


##1.8

trimesh_pbmc_umap <- ggplot(df_bin_centroids_pbmc, aes(x = x, y = y)) +
  geom_segment(data = tr_from_to_df_pbmc, aes(x = x_from, y = y_from, xend = x_to, yend = y_to)) +
  geom_point(size = 2, colour = "#33a02c") +
  coord_equal()

# ggplot(df_bin_centroids, aes(x = x, y = y)) +
# geom_point(size = 1, colour = "#33a02c") +
# geom_trimesh() +
# coord_equal()

trimesh_pbmc_umap <- trimesh_pbmc_umap +
  #ggtitle("(a)") +
  xlab(expression(C[x]^{(2)})) + ylab(expression(C[y]^{(2)})) +
  theme_light() +
  theme(legend.position = "none", plot.title = element_text(size = 5, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()#change legend key width
  )
# theme(axis.text = element_text(size = 5),
#       axis.title = element_text(size = 7))

trimesh_gr_pbmc_umap <- colour_long_edges(.data = distance_pbmc, benchmark_value = benchmark_pbmc,
                                          triangular_object = tr1_object_pbmc, distance_col = distance)

trimesh_gr_pbmc_umap <- trimesh_gr_pbmc_umap +
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
  )

trimesh_removed_pbmc_umap <- remove_long_edges(.data = distance_pbmc, benchmark_value = benchmark_pbmc,
                                               triangular_object = tr1_object_pbmc, distance_col = distance)

trimesh_removed_pbmc_umap <- trimesh_removed_pbmc_umap +
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
  )

############ Compute absolute residuals

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

library(ggbeeswarm)
error_df$group <- "1"
ggplot(error_df, aes(x = group, y = total)) +
  geom_quasirandom()+
  ylim(0, max(unlist(error_df$total))+ 0.5) + coord_flip()

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


#df$type <- as.factor(df$type)

#levels(df$type) <- c("Memory \nCD4 T", "Naive CD4 T", "CD14+ Mono",  "B", "CD8 T", "FCGR3A+ \n Mono", "NK", "M-MDSC\n-like", "CD27-CD+ \n Memory T", "DC")

df_b <- df_bin_pbmc |>
  dplyr::filter(hb_id %in% df_bin_centroids_pbmc$hexID) |>
  dplyr::select(-hb_id) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

df_exe <- dplyr::bind_rows(df_b, df)

distance_df_small_edges <- distance_pbmc %>%
  dplyr::filter(distance < benchmark_pbmc)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to, group = factor(df_exe$type , levels = c("error 2 or less", "error 2-10", "error 10-15", "error 15-20", "error 20-25", "error 25-30", "error 30-35", "error greter than 35", "model")
                         ), pointSize = 3,
                         levelColors = c("#b15928", "#1f78b4", "#ccebc5",
                                         "#fb9a99", "#cab2d6", "#ff7f00",
                                         "#ffed6f", "#e31a1c", "#33a02c"))


tour1_pbmc_umap <- show_langevitour(df_all_pbmc, df_bin_pbmc,
                                    df_bin_centroids_pbmc, benchmark_value = benchmark_pbmc,
                                    distance = distance_pbmc, distance_col = distance, col_start = "PC")

