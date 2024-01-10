library(readr)
library(umap)
library(dplyr)
library(rsample)

set.seed(20240110)

source("nldr_code.R", local = TRUE)

data <- read_csv(paste0(here::here(), "/data/s_curve.csv"))

data <- data |>
  mutate(ID = row_number())

data_split <- initial_split(data)
training_data <- training(data_split) |>
  arrange(ID)
test_data <- testing(data_split) |>
  arrange(ID)

UMAP_fit <- umap(training_data |> dplyr::select(-ID), n_neighbors = 50, n_components =  2)

UMAP_data <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data)[1:(ncol(UMAP_data))] <- paste0(rep("UMAP",(ncol(UMAP_data))), 1:(ncol(UMAP_data)))

UMAP_data <- UMAP_data |>
  mutate(ID = training_data$ID)

## Run only once
write_rds(UMAP_data, file = "data/s_curve/s_curve_umap.rds")

predict_UMAP_df <- predict(UMAP_fit, test_data |> dplyr::select(-ID)) |>
  as.data.frame()

names(predict_UMAP_df)[1:(ncol(predict_UMAP_df))] <- paste0(rep("UMAP",(ncol(predict_UMAP_df))), 1:(ncol(predict_UMAP_df)))

predict_UMAP_df <- predict_UMAP_df |>
  mutate(ID = test_data$ID)

## Run only once
write_rds(UMAP_data, file = "data/s_curve/s_curve_umap_predict.rds")

