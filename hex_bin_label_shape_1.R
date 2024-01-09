## Shape value issue
cell_area <- 0.5

num_bins_x <- calculate_effective_x_bins(.data = UMAP_data, x = UMAP1,
                                         cell_area = 0.5)

shape_val <- calculate_effective_shape_value(.data = UMAP_data,
                                             x = UMAP1, y = UMAP2)


#shape_val <- 1

hb_object <- hexbin::hexbin(x = UMAP_data |> pull(UMAP1),
                            y = UMAP_data |> pull(UMAP2),
                            xbins = num_bins_x, IDs = TRUE,
                            shape = shape_val)

hexdf_data <- tibble::tibble(tibble::as_tibble(hexbin::hcell2xy(hb_object)),  hexID = hb_object@cell, counts = hb_object@count/max(hb_object@count))

hex_grid <- expand.grid(UMAP_data |> pull(UMAP1), UMAP_data |> pull(UMAP2))

hb <- hexbin::hexbin(x = hex_grid |> pull(Var1),
                     y = hex_grid |> pull(Var2),
                     xbins = num_bins_x, IDs = TRUE,
                     shape = shape_val)

hexdf <- tibble::tibble(tibble::as_tibble(hexbin::hcell2xy(hb)),  hexID = hb@cell, counts = hb@count/max(hb@count))

hexdf <- hexdf %>%
  mutate(counts = ifelse((hexID %in% (hexdf_data |> pull(hexID))),counts, NA))

embedding_hb <- create_hexbin(.data = training_data, nldr_df = UMAP_data,
                              embedding_1 = UMAP1, embedding_2 = UMAP2,
                              num_bins = num_bins_x, shape_val = shape_val, apply_pca = FALSE)$df_new

full_grid <- full_join(hexdf, embedding_hb, by = c("hexID" = "hb_id"))



ggplot(full_grid, aes(x = x, y = y, fill = counts, hexID = hexID, label = hexID)) +
  geom_hex(stat = "identity", color = "#969696") +
  geom_label() +
  geom_point(aes(x = umap1, y = umap2), na.rm = TRUE) +
  scale_fill_viridis_c(na.value = "#ffffff") +
  ggtitle(paste0(" A = ", cell_area, ", b = (", hb_object@dimen[2], ", ", hb_object@dimen[1], "), ", "s = ", round(shape_val, 3))) +
  #coord_equal() +
  theme_light() +
  theme(legend.position="bottom", legend.direction="horizontal", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm')) +
  guides(fill = guide_colourbar(title = "Standardized count"))


## Full grid issue




sample_size <- 1500
cluster_size <- sample_size/5
df1 <- tibble::tibble(x1=rnorm(cluster_size, mean = 0, sd = 0.05), x2=rnorm(cluster_size, mean = 0, sd = 0.05), x3=rnorm(cluster_size, mean = 0, sd = 0.05), x4=rnorm(cluster_size, mean = 0, sd = 0.05), x5=rnorm(cluster_size, mean = 0, sd = 0.05), x6=rnorm(cluster_size, mean = 0, sd = 0.05), x7=rnorm(cluster_size, mean = 0, sd = 0.05), x8=rnorm(cluster_size, mean = 0, sd = 0.05), x9=rnorm(cluster_size, mean = 0, sd = 0.05), x10=rnorm(cluster_size, mean = 0, sd = 0.05))

df2 <- tibble::tibble(x1=rnorm(cluster_size, mean = 1, sd = 0.05), x2=rnorm(cluster_size, mean = 0, sd = 0.05), x3=rnorm(cluster_size, mean = 0, sd = 0.05), x4=rnorm(cluster_size, mean = 0, sd = 0.05), x5=rnorm(cluster_size, mean = 0, sd = 0.05), x6=rnorm(cluster_size, mean = 0, sd = 0.05), x7=rnorm(cluster_size, mean = 0, sd = 0.05), x8=rnorm(cluster_size, mean = 0, sd = 0.05), x9=rnorm(cluster_size, mean = 0, sd = 0.05), x10=rnorm(cluster_size, mean = 0, sd = 0.05))

df3 <- tibble::tibble(x1=rnorm(cluster_size, mean = 0, sd = 0.05), x2=rnorm(cluster_size, mean = 1, sd = 0.05), x3=rnorm(cluster_size, mean = 0, sd = 0.05), x4=rnorm(cluster_size, mean = 0, sd = 0.05), x5=rnorm(cluster_size, mean = 0, sd = 0.05), x6=rnorm(cluster_size, mean = 0, sd = 0.05), x7=rnorm(cluster_size, mean = 0, sd = 0.05), x8=rnorm(cluster_size, mean = 0, sd = 0.05), x9=rnorm(cluster_size, mean = 0, sd = 0.05), x10=rnorm(cluster_size, mean = 0, sd = 0.05))

df4 <- tibble::tibble(x1=rnorm(cluster_size, mean = 0, sd = 0.05), x2=rnorm(cluster_size, mean = 0, sd = 0.05), x3=rnorm(cluster_size, mean = 1, sd = 0.05), x4=rnorm(cluster_size, mean = 0, sd = 0.05), x5=rnorm(cluster_size, mean = 0, sd = 0.05), x6=rnorm(cluster_size, mean = 0, sd = 0.05), x7=rnorm(cluster_size, mean = 0, sd = 0.05), x8=rnorm(cluster_size, mean = 0, sd = 0.05), x9=rnorm(cluster_size, mean = 0, sd = 0.05), x10=rnorm(cluster_size, mean = 0, sd = 0.05))

df5 <- tibble::tibble(x1=rnorm(cluster_size, mean = 0, sd = 0.05), x2=rnorm(cluster_size, mean = 0, sd = 0.05), x3=rnorm(cluster_size, mean = 0, sd = 0.05), x4=rnorm(cluster_size, mean = 1, sd = 0.05), x5=rnorm(cluster_size, mean = 0, sd = 0.05), x6=rnorm(cluster_size, mean = 0, sd = 0.05), x7=rnorm(cluster_size, mean = 0, sd = 0.05), x8=rnorm(cluster_size, mean = 0, sd = 0.05), x9=rnorm(cluster_size, mean = 0, sd = 0.05), x10=rnorm(cluster_size, mean = 0, sd = 0.05))

df_2 <- dplyr::bind_rows(df1, df2, df3, df4, df5)




data_split <- initial_split(df_2)
training_data <- training(data_split)
test_data <- testing(data_split)



tSNE_data <- Fit_tSNE(training_data, with_seed = 20230531)



tSNE_example1 <- plot_tSNE_2D(tSNE_data) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())




cell_area <- 1

num_bins <- calculate_effective_x_bins(.data = tSNE_data, x = tSNE1,
                                         cell_area = 1)
shape_val <- 1



hb <- hexbin::hexbin(tSNE_data %>% pull(tSNE1), tSNE_data %>% pull(tSNE2), num_bins, IDs = TRUE, shape = shape_val)



hexbin <- create_hexbin(training_data, tSNE_data, tSNE1, tSNE2, num_bins = num_bins, shape_val = shape_val, apply_pca = FALSE)

df_all <- hexbin$df
hb_object <- hexbin$hb




hb_data <- hexbin::hexbin(x = tSNE_data |> pull(tSNE1),
                          y = tSNE_data |> pull(tSNE2),
                          xbins = num_bins, IDs = TRUE,
                          shape = shape_val)

hexdf_data <- tibble::tibble(tibble::as_tibble(hexbin::hcell2xy(hb_data)),  hexID = hb_data@cell, counts = hb_data@count/max(hb_data@count))

hex_grid <- expand.grid(tSNE_data |> pull(tSNE1), tSNE_data |> pull(tSNE2))

hex_grid_all <- expand.grid(min(tSNE_data |> pull(tSNE1)): max(tSNE_data |> pull(tSNE1)), min(tSNE_data |> pull(tSNE2)): max(tSNE_data |> pull(tSNE2)))
#hex_grid_all <- expand.grid(min(hex_grid |> pull(Var1)): max(hex_grid |> pull(Var1)), min(hex_grid |> pull(Var2)): max(hex_grid |> pull(Var2)))

hex_grid <- bind_rows(hex_grid, hex_grid_all) %>%
  distinct()

hb <- hexbin::hexbin(x = hex_grid |> pull(Var1),
                     y = hex_grid |> pull(Var2),
                     xbins = num_bins, IDs = TRUE,
                     shape = shape_val)

hexdf <- tibble::tibble(tibble::as_tibble(hexbin::hcell2xy(hb)),  hexID = hb@cell, counts = hb@count/max(hb@count))

hexdf <- hexdf %>%
  mutate(counts = ifelse((hexID %in% (hexdf_data |> pull(hexID))),counts, NA))

embedding_hb <- create_hexbin(.data = training_data, nldr_df = tSNE_data,
                              embedding_1 = tSNE1, embedding_2 = tSNE2,
                              num_bins = num_bins, shape_val = shape_val, apply_pca = FALSE)$df_new

full_grid <- full_join(hexdf, embedding_hb, by = c("hexID" = "hb_id"))

ggplot(full_grid, aes(x = x, y = y, fill = counts, hexID = hexID)) +
  geom_hex(stat = "identity", color = "#969696") +
  geom_point(aes(x = tsne1, y = tsne2), na.rm = TRUE) +
  scale_fill_viridis_c(na.value = "#ffffff") +
  ggtitle(paste0("A = ", 1 , ", b = (", hb_data@dimen[2], ", ", (hb_data@dimen[1] - 1), ")")) +
  theme_light() +
  #coord_equal() +
  theme(legend.position = "bottom", plot.title = element_text(size = 5, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm')) +
  guides(fill = guide_colourbar(title = "Standardized count"))



