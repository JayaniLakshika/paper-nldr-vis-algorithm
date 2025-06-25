library(dplyr)
library(snedata)
library(ggflowchart)
library(purrr) ## map function
library(gridExtra) ## for grid.arrange
library(rsample)
library(DT)
library(ggbeeswarm)
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

sample_size <- 2000
cluster_size <- sample_size/2
df1 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.1),
                      y=rnorm(cluster_size, mean = 0, sd = 0.1),
                      z=rnorm(cluster_size, mean = 0, sd = 0.1),
                      w=rnorm(cluster_size, mean = 0, sd = 0.1))

df2 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.1),
                      y=rnorm(cluster_size, mean = 0, sd = 0.1),
                      z=rnorm(cluster_size, mean = 0, sd = 0.1),
                      w=rnorm(cluster_size, mean = 1, sd = 0.1))

df_2 <- bind_rows(df1, df2)
df_2 <- df_2 |>
  rename(x1 = x, x2 = y, x3 = z, x4 = w)

df_2 <- df_2 |>
  mutate(ID = row_number())

write_rds(df_2, file = "data/one_gau_cluster/data_one_gau.rds")


data_split_sp <- initial_split(df_2)
training_data_5 <- training(data_split_sp) |>
  arrange(ID)
test_data_5 <- testing(data_split_sp) |>
  arrange(ID)

write_rds(training_data_5, file = "data/one_gau_cluster/data_one_gau_training.rds")
write_rds(test_data_5, file = "data/one_gau_cluster/data_one_gau_test.rds")

### tSNE
tSNE_data_gau <- Fit_tSNE(training_data_5 |> dplyr::select(-ID), opt_perplexity = calculate_effective_perplexity(training_data_5), with_seed = 20240110)

tSNE_data_gau <- tSNE_data_gau |>
  select(-ID) |>
  mutate(ID = training_data_5$ID)

plot_tSNE_2D(tSNE_data_gau)

write_rds(tSNE_data_gau, file = paste0("data/one_gau_cluster/tsne_data_one_gau_", calculate_effective_perplexity(training_data_5),".rds"))


### UMAP

UMAP_fit <- umap::umap(training_data_5 |> dplyr::select(-ID), n_neighbors = 15, n_components =  2)

UMAP_data_gau <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data_gau)[1:(ncol(UMAP_data_gau))] <- paste0(rep("UMAP",(ncol(UMAP_data_gau))), 1:(ncol(UMAP_data_gau)))

UMAP_data_gau <- UMAP_data_gau |>
  mutate(ID = training_data_5$ID)

write_rds(UMAP_data_gau, file = "data/one_gau_cluster/umap_data_one_gau.rds")

## predict umap embeddings

predict_UMAP_df <- predict(UMAP_fit, test_data_5 |> dplyr::select(-ID)) |>
  as.data.frame()

names(predict_UMAP_df)[1:(ncol(predict_UMAP_df))] <- paste0(rep("UMAP",(ncol(predict_UMAP_df))), 1:(ncol(predict_UMAP_df)))

predict_UMAP_df <- predict_UMAP_df |>
  mutate(ID = test_data_5$ID)

plot_UMAP_2D(UMAP_data_gau)

write_rds(predict_UMAP_df, file = "data/one_gau_cluster/umap_data_one_gau_predict.rds")


### TriMAP

tem_dir <- tempdir()

Fit_TriMAP_data(training_data_5 |> dplyr::select(-ID), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_TriMAP_values.csv")

Fit_TriMAP(as.integer(2), as.integer(5), as.integer(4), as.integer(3), path, path2)

TriMAP_data <- read_csv(path2)
TriMAP_data <- TriMAP_data |>
  mutate(ID = training_data_5$ID)
write_rds(TriMAP_data, file = "data/one_gau_cluster/trimap_data_one_gau.rds")


### PaCMAP

tem_dir <- tempdir()

Fit_PacMAP_data(training_data_5 |> dplyr::select(-ID), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_PaCMAP_values.csv")

Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.9, as.integer(2), path, path2)

PacMAP_data <- read_csv(path2)
PacMAP_data <- PacMAP_data |>
  mutate(ID = training_data_5$ID)

write_rds(PacMAP_data, file = "data/one_gau_cluster/pacmap_data_one_gau.rds")



### Phate

PHATE_data <- Fit_PHATE(training_data_5 |> dplyr::select(-ID), knn = 5, with_seed = 20240110)
PHATE_data <- PHATE_data |>
  select(PHATE1, PHATE2)
PHATE_data <- PHATE_data |>
  mutate(ID = training_data_5$ID)

write_rds(PHATE_data, file = "data/one_gau_cluster/phate_data_one_gau.rds")

