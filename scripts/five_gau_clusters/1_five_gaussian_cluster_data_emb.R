library(dplyr)
library(purrr) ## map function
library(rsample)
library(ggplot2)
library(readr)

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

sample_size <- 5000
cluster_size <- sample_size/5
df1 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.05), y=rnorm(cluster_size, mean = 0, sd = 0.05), z=rnorm(cluster_size, mean = 0, sd = 0.05), w=rnorm(cluster_size, mean = 0, sd = 0.05))
df2 <- tibble::tibble(x=rnorm(cluster_size, mean = 1, sd = 0.05), y=rnorm(cluster_size, mean = 0, sd = 0.05), z=rnorm(cluster_size, mean = 0, sd = 0.05), w=rnorm(cluster_size, mean = 0, sd = 0.05))
df3 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.05), y=rnorm(cluster_size, mean = 1, sd = 0.05), z=rnorm(cluster_size, mean = 0, sd = 0.05), w=rnorm(cluster_size, mean = 0, sd = 0.05))
df4 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.05), y=rnorm(cluster_size, mean = 0, sd = 0.05), z=rnorm(cluster_size, mean = 1, sd = 0.05), w=rnorm(cluster_size, mean = 0, sd = 0.05))
df5 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.05), y=rnorm(cluster_size, mean = 0, sd = 0.05), z=rnorm(cluster_size, mean = 0, sd = 0.05), w=rnorm(cluster_size, mean = 1, sd = 0.05))

df_2 <- bind_rows(df1, df2, df3, df4, df5)
df_2 <- df_2 |>
  rename(x1 = x, x2 = y, x3 = z, x4 = w)

df_2 <- df_2 |>
  mutate(ID = row_number())

write_rds(df_2, file = "data/five_gau_clusters/data_five_gau.rds")

## With cluster labels

df1 <- df1 |>
  mutate(cluster = "cluster1")
df2 <- df2 |>
  mutate(cluster = "cluster2")
df3 <- df3 |>
  mutate(cluster = "cluster3")
df4 <- df4 |>
  mutate(cluster = "cluster4")
df5 <- df5 |>
  mutate(cluster = "cluster5")

df_2n <- bind_rows(df1, df2, df3, df4, df5)
df_2n <- df_2n |>
  rename(x1 = x, x2 = y, x3 = z, x4 = w)

df_2n <- df_2n |>
  mutate(ID = row_number())

write_rds(df_2n, file = "data/five_gau_clusters/data_five_gau_with_labels.rds")


# data_split_sp <- initial_split(df_2)
# training_data_5 <- training(data_split_sp) |>
#   arrange(ID)
# test_data_5 <- testing(data_split_sp) |>
#   arrange(ID)
#
# write_rds(training_data_5, file = "data/five_gau_clusters/data_five_gau_training.rds")
# write_rds(test_data_5, file = "data/five_gau_clusters/data_five_gau_test.rds")

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

## predict umap embeddings

predict_UMAP_df <- predict(UMAP_fit, test_data_5 |> dplyr::select(-ID)) |>
  as.data.frame()

names(predict_UMAP_df)[1:(ncol(predict_UMAP_df))] <- paste0(rep("UMAP",(ncol(predict_UMAP_df))), 1:(ncol(predict_UMAP_df)))

predict_UMAP_df <- predict_UMAP_df |>
  mutate(ID = test_data_5$ID)

plot_UMAP_2D(UMAP_data_gau)

write_rds(predict_UMAP_df, file = "data/five_gau_clusters/umap_data_five_gau_predict.rds")


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

