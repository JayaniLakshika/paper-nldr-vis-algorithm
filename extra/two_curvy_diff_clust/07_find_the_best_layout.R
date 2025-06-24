library(ggplot2)
library(dplyr)
library(readr)

conflicted::conflicts_prefer(dplyr::filter)

error_two_c_shaped_clust_umap <- read_rds("data/two_c_shaped_clust/error_two_c_shaped_clust_umap.rds")
error_two_c_shaped_clust_tsne <- read_rds("data/two_c_shaped_clust/error_two_c_shaped_clust_tsne.rds")
error_two_c_shaped_clust_phate <- read_rds("data/two_c_shaped_clust/error_two_c_shaped_clust_phate.rds")
error_two_c_shaped_clust_trimap <- read_rds("data/two_c_shaped_clust/error_two_c_shaped_clust_trimap.rds")
error_two_c_shaped_clust_pacmap <- read_rds("data/two_c_shaped_clust/error_two_c_shaped_clust_pacmap.rds")

error_two_c_shaped_clust <- bind_rows(error_two_c_shaped_clust_umap,
                                   error_two_c_shaped_clust_tsne,
                                   error_two_c_shaped_clust_phate,
                                   error_two_c_shaped_clust_trimap,
                                   error_two_c_shaped_clust_pacmap)

error_two_c_shaped_clust <- error_two_c_shaped_clust |>
  mutate(a1 = round(a1, 2)) |>
  filter(bin1 >= 5) |>
  group_by(method, a1) |>
  filter(MSE == min(MSE)) |>
  ungroup()

error_plot_two_c_shaped_clust <- ggplot(error_two_c_shaped_clust,
                                     aes(x = a1,
                                         y = MSE,
                                         colour = method)) +
  geom_point(size = 0.8) +
  geom_line(linewidth = 0.3) +
  # geom_vline(xintercept = 15, linetype="solid",
  #            color = "black", linewidth=0.8, alpha = 0.5) +
  scale_x_continuous(breaks = sort(unique(error_two_c_shaped_clust$a1))[seq(1, length(unique(error_two_c_shaped_clust$a1)), by = 5)]) +
  scale_color_manual(values=c('#e41a1c','#ff7f00','#4daf4a', "#a65628",'#636363')) +
  scale_y_log10() +
  ylab("log(MSE)") +
  xlab(expression(paste("binwidth (", a[1], ")"))) +
  theme_minimal() +
  theme(panel.border = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 12, hjust = 0.5, vjust = -0.5),
        axis.ticks.x = element_line(),
        axis.ticks.y = element_line(),
        legend.position = "bottom",
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7))

error_plot_two_c_shaped_clust

## PaCMAP is best because in tSNE the non-linear curve is separated.
