library(dplyr)
library(readr)
library(ggplot2)

tsne_df <- read_rds(here::here("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_30.rds"))
tsne_df2 <- read_rds(here::here("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_62.rds"))
umap_df <- read_rds(here::here("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_umap_n-neigbors_15_min-dist_0.1.rds"))
phate_df <- read_rds(here::here("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_phate_knn_5.rds"))
trimap_df <- read_rds(here::here("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_trimap_n-inliers_12_n-outliers_4_n-random_3.rds"))
pacmap_df <- read_rds(here::here("data/two_non_linear_diff_shaped_close_clusters/two_non_linear_diff_shaped_close_clusters_pacmap_n-neighbors_10_init_random_MN-ratio_0.5_FP-ratio_2.rds"))


tsne_df <- tsne_df |>
  mutate(method = "tSNE")
names(tsne_df) <- c("emb_1", "emb_2", "method")

tsne_df2 <- tsne_df2 |>
  mutate(method = "tSNE2")
names(tsne_df2) <- c("emb_1", "emb_2", "method")

umap_df <- umap_df |>
  mutate(method = "UMAP")
names(umap_df) <- c("emb_1", "emb_2", "method")

phate_df <- phate_df |>
  mutate(method = "PHATE")
names(phate_df) <- c("emb_1", "emb_2", "method")

trimap_df <- trimap_df |>
  mutate(method = "TriMAP")
names(trimap_df) <- c("emb_1", "emb_2", "method")

pacmap_df <- pacmap_df |>
  mutate(method = "PaCMAP")
names(pacmap_df) <- c("emb_1", "emb_2", "method")

nldr_df <- bind_rows(tsne_df,
                     tsne_df2,
                     umap_df,
                     phate_df,
                     pacmap_df,
                     trimap_df)


ggplot(data = nldr_df,
       aes(
         x = emb_1,
         y = emb_2
       )) +
  geom_point(
    alpha = 0.5
  ) +
  facet_wrap(~factor(method, levels = c("tSNE", "tSNE2", "UMAP", "PHATE", "TriMAP", "PaCMAP")),
             scales = "free")
