library(dplyr)
library(purrr) ## map function
library(rsample)
library(ggplot2)
library(readr)
library(geozoo)
library(mvtnorm)

library(Rtsne)
library(umap)
library(phateR)
library(reticulate)
library(patchwork)

library(grid)


set.seed(20240110)

use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

reticulate::source_python(here::here("scripts/function_scripts/Fit_PacMAP_code.py"))
reticulate::source_python(here::here("scripts/function_scripts/Fit_TriMAP_code.py"))

source(here::here("scripts/nldr_code.R"))

p2 <- 4
vert2 <- simplex(p2)$points

df <- rmvnorm(5000, mean = rep(0, 4), sigma = matrix(c(c(0.1, 0.5, 1, 0.3),
                                                       c(1, 0.05, 0.8, 0.2),
                                                       c(1, 0.05, 0.8, 2),
                                                       c(1.5, 0.05, 0.8, 0.2)), 4, 4))/6 #diag(4) * 0.05
df[1:1000,] <- df[1:1000,] + matrix(rep(vert2[1,], 1000), ncol=4, byrow=T)
df[1001:2000,] <- df[1001:2000,] + matrix(rep(vert2[2,], 1000), ncol=4, byrow=T)
df[2001:3000,] <- df[2001:3000,] + matrix(rep(vert2[3,], 1000), ncol=4, byrow=T)
df[3001:4000,] <- df[3001:4000,] + matrix(rep(vert2[4,], 1000), ncol=4, byrow=T)
df[4001:5000,] <- df[4001:5000,] + matrix(rep(vert2[5,], 1000), ncol=4, byrow=T)
df_2 <- as_tibble(df)

df_2 <- df_2 |>
  mutate(across(everything(), ~ (. - mean(.)) / sd(.)))

names(df_2) <- paste0("x", 1:4)

df_2 <- df_2 |>
  mutate(ID = row_number())

write_rds(df_2, file = "data/five_gau_clusters/data_five_gau.rds")

### tSNE
tSNE_data_gau <- Fit_tSNE(df_2 |> dplyr::select(-ID), opt_perplexity = calculate_effective_perplexity(df_2), with_seed = 20240110)

tSNE_data_gau <- tSNE_data_gau |>
  select(-ID) |>
  mutate(ID = df_2$ID)

plot_tSNE_2D(tSNE_data_gau)

write_rds(tSNE_data_gau, file = paste0("data/five_gau_clusters/tsne_data_five_gau_", calculate_effective_perplexity(df_2),".rds"))
#write_rds(tSNE_data_gau, file = paste0("data/five_gau_clusters/tsne_data_five_gau_30.rds"))


### UMAP

UMAP_fit <- umap(df_2 |> dplyr::select(-ID), n_neighbors = 15, n_components =  2)

UMAP_data_gau <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data_gau)[1:(ncol(UMAP_data_gau))] <- paste0(rep("UMAP",(ncol(UMAP_data_gau))), 1:(ncol(UMAP_data_gau)))

UMAP_data_gau <- UMAP_data_gau |>
  mutate(ID = df_2$ID)

write_rds(UMAP_data_gau, file = "data/five_gau_clusters/umap_data_five_gau.rds")

### TriMAP

tem_dir <- tempdir()

Fit_TriMAP_data(df_2 |> dplyr::select(-ID), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_TriMAP_values.csv")

#Fit_TriMAP(as.integer(2), as.integer(5), as.integer(4), as.integer(3), path, path2)
Fit_TriMAP(as.integer(2), as.integer(12), as.integer(4), as.integer(3), path, path2)

TriMAP_data <- read_csv(path2)
TriMAP_data <- TriMAP_data |>
  mutate(ID = df_2$ID)
write_rds(TriMAP_data, file = "data/five_gau_clusters/trimap_data_five_gau.rds")


### PaCMAP

tem_dir <- tempdir()

Fit_PacMAP_data(df_2 |> dplyr::select(-ID), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_PaCMAP_values.csv")

# Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.9, as.integer(2), path, path2)
Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.5, as.integer(2), path, path2)

PacMAP_data <- read_csv(path2)
PacMAP_data <- PacMAP_data |>
  mutate(ID = df_2$ID)

write_rds(PacMAP_data, file = "data/five_gau_clusters/pacmap_data_five_gau.rds")



### Phate

PHATE_data <- Fit_PHATE(df_2 |> dplyr::select(-ID), knn = 5, with_seed = 20240110)
PHATE_data <- PHATE_data |>
  select(PHATE1, PHATE2)
PHATE_data <- PHATE_data |>
  mutate(ID = df_2$ID)

write_rds(PHATE_data, file = "data/five_gau_clusters/phate_data_five_gau.rds")

