library(readr)
library(umap)
library(dplyr)


training_data <- read_rds("data/s_curve/s_curve_training.rds")
test_data <- read_rds("data/s_curve/s_curve_test.rds")

## UMAP with n_neighbors 7

UMAP_fit <- umap(training_data |> dplyr::select(-ID), n_neighbors = 7, n_components =  2)

UMAP_data <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data)[1:(ncol(UMAP_data))] <- paste0(rep("UMAP",(ncol(UMAP_data))), 1:(ncol(UMAP_data)))

UMAP_data <- UMAP_data |>
  mutate(ID = training_data$ID)

## Run only once
write_rds(UMAP_data, file = "data/s_curve/s_curve_umap_7.rds")

predict_UMAP_df <- predict(UMAP_fit, test_data |> dplyr::select(-ID)) |>
  as.data.frame()

names(predict_UMAP_df)[1:(ncol(predict_UMAP_df))] <- paste0(rep("UMAP",(ncol(predict_UMAP_df))), 1:(ncol(predict_UMAP_df)))

predict_UMAP_df <- predict_UMAP_df |>
  mutate(ID = test_data$ID)

## Run only once
write_rds(UMAP_data, file = "data/s_curve/s_curve_umap_7_predict.rds")




## UMAP with n_neighbors 15

UMAP_fit <- umap(training_data |> dplyr::select(-ID), n_neighbors = 15, n_components =  2)

UMAP_data <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data)[1:(ncol(UMAP_data))] <- paste0(rep("UMAP",(ncol(UMAP_data))), 1:(ncol(UMAP_data)))

UMAP_data <- UMAP_data |>
  mutate(ID = training_data$ID)

## Run only once
write_rds(UMAP_data, file = "data/s_curve/s_curve_umap_15.rds")

predict_UMAP_df <- predict(UMAP_fit, test_data |> dplyr::select(-ID)) |>
  as.data.frame()

names(predict_UMAP_df)[1:(ncol(predict_UMAP_df))] <- paste0(rep("UMAP",(ncol(predict_UMAP_df))), 1:(ncol(predict_UMAP_df)))

predict_UMAP_df <- predict_UMAP_df |>
  mutate(ID = test_data$ID)

## Run only once
write_rds(UMAP_data, file = "data/s_curve/s_curve_umap_15_predict.rds")


## UMAP with n_neighbors 32

UMAP_fit <- umap(training_data |> dplyr::select(-ID), n_neighbors = 32, n_components =  2)

UMAP_data <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data)[1:(ncol(UMAP_data))] <- paste0(rep("UMAP",(ncol(UMAP_data))), 1:(ncol(UMAP_data)))

UMAP_data <- UMAP_data |>
  mutate(ID = training_data$ID)

## Run only once
write_rds(UMAP_data, file = "data/s_curve/s_curve_umap_32.rds")

predict_UMAP_df <- predict(UMAP_fit, test_data |> dplyr::select(-ID)) |>
  as.data.frame()

names(predict_UMAP_df)[1:(ncol(predict_UMAP_df))] <- paste0(rep("UMAP",(ncol(predict_UMAP_df))), 1:(ncol(predict_UMAP_df)))

predict_UMAP_df <- predict_UMAP_df |>
  mutate(ID = test_data$ID)

## Run only once
write_rds(UMAP_data, file = "data/s_curve/s_curve_umap_32_predict.rds")

