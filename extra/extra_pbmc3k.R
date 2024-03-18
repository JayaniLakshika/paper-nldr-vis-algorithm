```{r}
#| warning: false
#| echo: false

## Import data
training_data_pbmc <- read_rds("data/pbmc/pbmc_3k_festem/pbmc_pca_50.rds")
training_data_pbmc <- training_data_pbmc[, 1:15] |>
  mutate(ID = 1:NROW(training_data_pbmc))

umap_pbmc <- read_rds("data/pbmc/pbmc_3k_festem/pbmc_umap.rds")
umap_pbmc_scaled <- as.data.frame(do.call(cbind, gen_scaled_data(data = umap_pbmc,
                                                                 x = "UMAP1", y = "UMAP2"))) |>
  dplyr::mutate(cell_label = umap_pbmc$cell_label) |>
  dplyr::rename(c("UMAP1" = "scaled_UMAP1",
                  "UMAP2" = "scaled_UMAP2")) |>
  dplyr::mutate(ID = 1:NROW(umap_pbmc))

class_avg <- umap_pbmc_scaled |>
  group_by(cell_label) |>
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2)
  ) |>
  mutate(cell_label = as.numeric(cell_label)) |>
  arrange(cell_label)

class_avg$cell_label <- as.factor(class_avg$cell_label)

levels(class_avg$cell_label) <- c("Memory \nCD4 T", "Naive CD4 T", "CD14+ Mono",  "B", "CD8 T", "FCGR3A+ \n Mono", "NK", "M-MDSC\n-like", "CD27-CD+ \n Memory T", "DC")

plot_list2_pbmc <- umap_pbmc_scaled |>
  ggplot(aes(x = UMAP1,
             y = UMAP2, color = cell_label))+
  geom_point(alpha=0.5) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + #ggtitle("(a)") +
  theme_linedraw() +
  theme(legend.position = "none", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm')) +
  scale_color_manual(values=c("#b15928", "#1f78b4", "#cab2d6", "#ccebc5", "#fb9a99", "#e31a1c", "#6a3d9a", "#ff7f00", "#ffed6f", "#fdbf6f", "#ffff99", "#a6cee3", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#b2df8a", "#bc80bd", "#33a02c", "#ccebc5", "#ffed6f", "#000000", "#bdbdbd")) +
  geom_text(aes(x=UMAP1,y = UMAP2,label = cell_label), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 2) +
  annotate(geom = 'text', label = 'a', x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 3)
```

```{r}
#| warning: false
#| echo: false

## Import data
training_data_pbmc <- read_rds("data/pbmc/pbmc_3k_festem/pbmc_pca_50.rds")
training_data_pbmc <- training_data_pbmc[, 1:15] |>
  mutate(ID = 1:NROW(training_data_pbmc))

umap_pbmc <- read_rds("data/pbmc/pbmc_umap_5_min_dist_0.99_metric_cosine.rds")
umap_pbmc_scaled <- as.data.frame(do.call(cbind, gen_scaled_data(data = umap_pbmc,
                                                                 x = "UMAP1", y = "UMAP2"))) |>
  dplyr::mutate(cell_label = umap_pbmc$cell_label) |>
  dplyr::rename(c("UMAP1" = "scaled_UMAP1",
                  "UMAP2" = "scaled_UMAP2")) |>
  dplyr::mutate(ID = 1:NROW(umap_pbmc))


class_avg <- umap_pbmc_scaled |>
  group_by(cell_label) |>
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2)
  ) |>
  mutate(cell_label = as.numeric(cell_label)) |>
  arrange(cell_label)

class_avg$cell_label <- as.factor(class_avg$cell_label)

levels(class_avg$cell_label) <- c("Memory \nCD4 T", "Naive CD4 T", "CD14+ Mono",  "B", "CD8 T", "FCGR3A+ \n Mono", "NK", "M-MDSC\n-like", "CD27-CD+ \n Memory T", "DC")

plot_list2_pbmc <- umap_pbmc_scaled |>
  ggplot(aes(x = UMAP1,
             y = UMAP2, color = cell_label))+
  geom_point(alpha=0.5) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + #ggtitle("(a)") +
  theme_linedraw() +
  theme(legend.position = "none", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm')) +
  scale_color_manual(values=c("#b15928", "#1f78b4", "#cab2d6", "#ccebc5", "#fb9a99", "#e31a1c", "#6a3d9a", "#ff7f00", "#ffed6f", "#fdbf6f", "#ffff99", "#a6cee3", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#b2df8a", "#bc80bd", "#33a02c", "#ccebc5", "#ffed6f", "#000000", "#bdbdbd")) +
  geom_text(aes(x=UMAP1,y = UMAP2,label = cell_label), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 2) +
  annotate(geom = 'text', label = 'a', x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5, size = 3)
```

