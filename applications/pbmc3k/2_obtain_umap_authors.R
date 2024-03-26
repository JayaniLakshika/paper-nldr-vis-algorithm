library(readr)
library(uwot)
library(dplyr)
set.seed(20240110)

## Select PCs
pbmc_pca <- read_rds("data/pbmc3k/pbmc_pca_50.rds")
pbmc_pca <- pbmc_pca[, 1:9] ## By looking at PCA scree plot in https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

## Obtain UMAP
umap_pbmc <- umap(pbmc_pca, n_neighbors = 30, n_components =  2,
                 metric = "cosine", min_dist = 0.3, init = "spca")

umap_pbmc <- umap_pbmc |>
  as.data.frame()

names(umap_pbmc)[1:(ncol(umap_pbmc))] <- paste0(rep("UMAP",(ncol(umap_pbmc))), 1:(ncol(umap_pbmc)))

umap_pbmc <- umap_pbmc |>
  mutate(ID = 1:NROW(umap_pbmc))

## Run only once
write_rds(umap_pbmc, file = "data/pbmc3k/pbmc_umap_30_min_dist_0.3.rds")
