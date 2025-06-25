library(dplyr)
library(readr)
library(ggplot2)

tsne_df <- read_rds(here::here("data/two_c_shaped_clust/two_c_shaped_clust_tsne_perplexity_30.rds"))
umap_df <- read_rds(here::here("data/two_c_shaped_clust/two_c_shaped_clust_umap_n-neigbors_15_min-dist_0.1.rds"))
phate_df <- read_rds(here::here("data/two_c_shaped_clust/two_c_shaped_clust_phate_knn_5.rds"))
trimap_df <- read_rds(here::here("data/two_c_shaped_clust/two_c_shaped_clust_trimap_n-inliers_12_n-outliers_4_n-random_3.rds"))
pacmap_df <- read_rds(here::here("data/two_c_shaped_clust/two_c_shaped_clust_pacmap_n-neighbors_10_init_random_MN-ratio_0.5_FP-ratio_2.rds"))


tsne_df <- tsne_df |>
  mutate(method = "tSNE")
names(tsne_df) <- c("emb_1", "emb_2", "method")

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
  facet_wrap(~factor(method, levels = c("tSNE", "UMAP", "PHATE", "TriMAP", "PaCMAP")),
             scales = "free")
