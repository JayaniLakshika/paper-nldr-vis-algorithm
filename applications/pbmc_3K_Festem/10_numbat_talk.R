## To import necessary packages
library(readr)
library(ggplot2)
library(langevitour)
library(detourr)

## Import UMAP data (Suggested by the author)
umap_pbmc <- read_rds("data/pbmc/pbmc_3k_festem/pbmc_umap.rds")

## Visualise UMAP
plot_list2_pbmc <- umap_pbmc |>
  ggplot(aes(x = UMAP1,
             y = UMAP2))+
  geom_point(alpha=0.5) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
  theme_linedraw() +
  theme(legend.position = "none", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'))

plot_list2_pbmc


## Import PBMC3k PCs
training_data_pbmc <- read_rds("data/pbmc/pbmc_3k_festem/pbmc_pca_50.rds")
## Visualise in high-D
langevitour(training_data_pbmc[, 1:15], levelColors = "black")

## Spin-brush approach with detour
detour(training_data_pbmc,
       tour_aes(projection = PC_1:PC_15)) |>
  tour_path(grand_tour(2), fps = 30,
            max_bases=100) |>
  show_scatter(alpha = 0.7,
               axes = FALSE)
## Save the clusters

## UMAP with my parameter choice
umap_pbmc <- read_rds("data/pbmc/pbmc_umap_5_min_dist_0.99_metric_cosine.rds")

## Visualise UMAP
plot_list2_pbmc <- umap_pbmc |>
  ggplot(aes(x = UMAP1,
             y = UMAP2))+
  geom_point(alpha=0.5) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
  theme_linedraw() +
  theme(legend.position = "none", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'))

plot_list2_pbmc

## Brush clusters
UMAP_pbmc_with_result1_label <- umap_pbmc |>
  dplyr::mutate(result1_class_label = case_when(
    ((UMAP1 >=-12) & (UMAP1 <= -5) & (UMAP2 >= 6) & (UMAP2 <= 18)) ~ "cluster 1",
    ((UMAP1 >-5) & (UMAP1 <= 7) & (UMAP2 >= 4.8) & (UMAP2 <= max(UMAP2))) ~ "cluster 2",
    ((UMAP1 >=6) & (UMAP1 <= 15) & (UMAP2 >= -4) & (UMAP2 < 4.8)) ~ "cluster 3",
    ((UMAP1 >=-3) & (UMAP1 <= 6.6) & (UMAP2 >= min(UMAP2)) & (UMAP2 <= 4.8)) ~ "cluster 4",
    ((UMAP1 >=min(UMAP1))& (UMAP1 < -3) & (UMAP2 >= -8) & (UMAP2 <= -1)) ~ "cluster 5",
    ((UMAP1 >=min(UMAP1)) & (UMAP1 <= -7) & (UMAP2 >= -3) & (UMAP2 <= 3)) ~ "cluster 6",
    .default = "no class"
  ))

### Reassign some labels

UMAP_pbmc_with_result1_label <- UMAP_pbmc_with_result1_label |>
  dplyr::mutate(result1_class_label = if_else((ID %in% c(2619, 1158)), "cluster 3", result1_class_label)) |>
  dplyr::mutate(result1_class_label = if_else((ID %in% c(948, 285, 158, 268, 731, 1152, 1372, 1651, 2010,
                                                         150, 2304, 2546, 487, 2176, 2157, 2451, 264, 63,
                                                         1234, 947, 1126, 1699, 1232, 1492, 1488, 2040, 2550,
                                                         1313, 361, 372, 1515, 1264, 1307, 1726, 517, 2568, 2274,
                                                         591, 1703, 1097, 430, 1556, 2133, 1232, 1194, 1888, 1796,
                                                         494, 1194)), "cluster 6", result1_class_label))


## Import detour results
detour_results8 <- read_csv("data/pbmc/pbmc_3k_festem/detourr_export8.csv")

### Rename labels

detour_results8 <- detour_results8 |>
  dplyr::mutate(result2_class_label = case_when(
    colour == "000000" ~ "cluster 4",
    colour == "49df20" ~ "cluster 2",
    colour == "619cff" ~ "cluster 3",
    colour == "d82cd3" ~ "cluster 1",
    colour == "df2020" ~ "cluster 5",
    colour == "ffbd61" ~ "cluster 6",
    .default = "no cluster"

  ))

df_res_1_2_label <- bind_cols(UMAP_pbmc_with_result1_label, detour_results8)

result1_2d_layout <- df_res_1_2_label |>
  ggplot(aes(x = UMAP1,
             y = UMAP2, color = result1_class_label))+
  geom_point(alpha=0.5) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
  theme_linedraw() +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm')) +
  scale_color_manual(values=c("#d82cd3", "#49df20", "#619cff", "#000000", "#df2020", "#ffbd61", "yellow")) +
  # geom_text(aes(x=UMAP1,y = UMAP2,label = cell_label), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 2) +
  annotate(geom = 'text', label = 'colored manually in 2D', x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3) +
  labs(colour = "cluster labels")

result2_2d_layout <- df_res_1_2_label |>
  ggplot(aes(x = UMAP1,
             y = UMAP2, color = result2_class_label))+
  geom_point(alpha=0.5) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
  theme_linedraw() +
  theme(legend.position = "bottom", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm')) +
  scale_color_manual(values=c("#d82cd3", "#49df20", "#619cff", "#000000", "#df2020", "#ffbd61", "yellow")) +
  # geom_text(aes(x=UMAP1,y = UMAP2,label = cell_label), data = class_avg,inherit.aes = F, color = "black",fontface = "bold",size = 2) +
  annotate(geom = 'text', label = 'colored in high-dimensions', x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3) +
  labs(colour = "cluster labels")

result1_2d_layout + result2_2d_layout +
  plot_layout(guides='collect', ncol=2) &
  theme(legend.position='bottom')

table(df_res_1_2_label$result1_class_label, df_res_1_2_label$result2_class_label)
