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
library(langevitour) ## zeiselPC is here

library(grid)

set.seed(20240110)

use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

reticulate::source_python(paste0(here::here(), "/examples/function_scripts/Fit_PacMAP_code.py"))
reticulate::source_python(paste0(here::here(), "/examples/function_scripts/Fit_TriMAP_code.py"))

source("nldr_code.R", local = TRUE)

data(zeiselPC)

zeiselPC <- zeiselPC[, 2:8]

zeiselPC <- zeiselPC |>
  mutate(ID = row_number())

write_rds(zeiselPC, "data/zeisel/zeisel.rds")

data_split <- initial_split(zeiselPC)
training_data <- training(data_split) |>
  arrange(ID)
test_data <- testing(data_split) |>
  arrange(ID)

write_rds(training_data, file = "data/zeisel/zeisel_training.rds")
write_rds(test_data, file = "data/zeisel/zeisel_test.rds")

## tSNE

tSNE_data <- Fit_tSNE(training_data |> select(-ID), opt_perplexity = 30, with_seed = 20240110)

tSNE_data <- tSNE_data |>
  select(-ID) |>
  mutate(ID = training_data$ID)

write_rds(tSNE_data, file = paste0("data/zeisel/zeisel_tsne_", 30, ".rds"))

## UMAP

UMAP_fit <- umap(training_data |> dplyr::select(-ID), n_neighbors = 15, n_components =  2)

UMAP_data <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data)[1:(ncol(UMAP_data))] <- paste0(rep("UMAP",(ncol(UMAP_data))), 1:(ncol(UMAP_data)))

UMAP_data <- UMAP_data |>
  mutate(ID = training_data$ID)

## Run only once
write_rds(UMAP_data, file = "data/zeisel/zeisel_umap.rds")

predict_UMAP_df <- predict(UMAP_fit, test_data |> dplyr::select(-ID)) |>
  as.data.frame()

names(predict_UMAP_df)[1:(ncol(predict_UMAP_df))] <- paste0(rep("UMAP",(ncol(predict_UMAP_df))), 1:(ncol(predict_UMAP_df)))

predict_UMAP_df <- predict_UMAP_df |>
  mutate(ID = test_data$ID)

## Run only once
write_rds(UMAP_data, file = "data/zeisel/zeisel_umap_predict.rds")

## PHATE

PHATE_data <- Fit_PHATE(training_data |> dplyr::select(-ID), knn = 5, with_seed = 20240110)
PHATE_data <- PHATE_data |>
  select(PHATE1, PHATE2)
PHATE_data <- PHATE_data |>
  mutate(ID = training_data$ID)

write_rds(PHATE_data, file = "data/zeisel/zeisel_phate.rds")


## TriMAP

tem_dir <- tempdir()

Fit_TriMAP_data(training_data |> dplyr::select(-ID), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_TriMAP_values.csv")

Fit_TriMAP(as.integer(2), as.integer(5), as.integer(4), as.integer(3), path, path2)

TriMAP_data <- read_csv(path2)
TriMAP_data <- TriMAP_data |>
  mutate(ID = training_data$ID)
write_rds(TriMAP_data, file = "data/zeisel/zeisel_trimap.rds")

## PaCMAP


tem_dir <- tempdir()

Fit_PacMAP_data(training_data |> dplyr::select(-ID), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_PaCMAP_values.csv")

Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.9, as.integer(2), path, path2)

PacMAP_data <- read_csv(path2)
PacMAP_data <- PacMAP_data |>
  mutate(ID = training_data$ID)

write_rds(PacMAP_data, file = "data/zeisel/zeisel_pacmap.rds")

