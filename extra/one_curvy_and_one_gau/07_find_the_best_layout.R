library(ggplot2)
library(dplyr)
library(readr)

conflicted::conflicts_prefer(dplyr::filter)

error_mnist_umap <- read_rds("data/one_curvy_one_gau_clust/error_one_curvy_one_gau_clust_tsne.rds")
error_mnist_tsne <- read_rds("data/one_curvy_one_gau_clust/error_one_curvy_one_gau_clust_umap.rds")
error_mnist_phate <- read_rds("data/one_curvy_one_gau_clust/error_one_curvy_one_gau_clust_phate.rds")
error_mnist_trimap <- read_rds("data/one_curvy_one_gau_clust/error_one_curvy_one_gau_clust_trimap.rds")
error_mnist_pacmap <- read_rds("data/one_curvy_one_gau_clust/error_one_curvy_one_gau_clust_pacmap.rds")

error_mnist <- bind_rows(error_mnist_umap,
                         error_mnist_tsne,
                         error_mnist_phate,
                         error_mnist_trimap,
                         error_mnist_pacmap)

error_mnist <- error_mnist |>
  mutate(a1 = round(a1, 2)) |>
  filter(bin1 >= 5) |>
  group_by(method, a1) |>
  filter(MSE == min(MSE)) |>
  ungroup()

error_plot_mnist <- ggplot(error_mnist,
                           aes(x = a1,
                               y = MSE,
                               colour = method)) +
  geom_point(size = 0.8) +
  geom_line(linewidth = 0.3) +
  # geom_vline(xintercept = 15, linetype="solid",
  #            color = "black", linewidth=0.8, alpha = 0.5) +
  scale_x_continuous(breaks = sort(unique(error_mnist$a1))[append(seq(1, length(unique(error_mnist$a1)), by = 5), 24)]) +
  scale_color_manual(values=c('#e41a1c','#377eb8','#4daf4a', "#ff7f00",'#a65628')) +
  scale_y_log10() +
  ylab("MSE") +
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

error_plot_mnist
