## This script is to compare the prediction approaches...

library(dplyr)
library(tibble)
library(rsample)
library(readr)
library(ggplot2)
library(conflicted)
library(quollr)

library(Rtsne)
library(umap) #predit for uwot not working
library(phateR)
library(reticulate)

set.seed(20240110)

conflicts_prefer(umap::umap)
conflicts_prefer(dplyr::filter)

use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

#reticulate::source_python(paste0(here::here("R/function_scripts/Fit_PacMAP_code.py")))
#reticulate::source_python(paste0(here::here("R/function_scripts/Fit_TriMAP_code.py")))

#source(here::here("R/nldr_code.R"), local = TRUE)

data <- read_rds(here::here("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_data.rds"))

################################################################################


# ##UMAP on the whole data (training + test)
#
# n_neighbors <- 15
# min_dist <- 0.1
#
# # Create a config list with the desired parameters
# umap_config <- umap.defaults
# umap_config$n_neighbors <- n_neighbors      # Set the number of neighbors
# umap_config$n_components <- 2    # Set the number of output dimensions (typically 2 or 3)
# umap_config$min_dist <- min_dist
#
# UMAP_fit <- umap(data, config = umap_config)
#
# UMAP_data <- UMAP_fit$layout |>
#   as_tibble()
#
# names(UMAP_data) <- c("UMAP1", "UMAP2")

## To split the data

## Data
UMAP_data <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_umap_n-neigbors_15_min-dist_0.1.rds")

UMAP_data <- UMAP_data |>
  mutate(ID = 1:NROW(UMAP_data)) |>
  mutate(cluster = if_else(ID <= 1000, "cluster1", "cluster2"))

split <- initial_split(UMAP_data, prop = 0.7, strata = cluster)
training_umap_two_curvy <- training(split)
test_umap_two_curvy <- testing(split)

### Predict the test set with UMAP model

# predict_UMAP_df <- predict(UMAP_fit, test_data_two_curvy[, 1:7]) |>
#   as_tibble()
#
# names(predict_UMAP_df) <- c("UMAP1", "UMAP2")

predict_UMAP_df <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_umap_predict_test.rds")


## To split data as well
data_two_curvy <- data |>
  mutate(ID = 1:NROW(data))

training_data_two_curvy <- data_two_curvy |>
  filter(ID %in% training_umap_two_curvy$ID)

test_data_two_curvy <- data_two_curvy |>
  filter(ID %in% test_umap_two_curvy$ID)

### quollr on training data using UMAP layout

## Fit the model

num_bins_x_two_curvy <- 54

algo_obj_two_curvy <- fit_highd_model(
  highd_data = training_data_two_curvy[, 1:8],
  nldr_data = training_umap_two_curvy[, 1:3],
  b1 = num_bins_x_two_curvy,
  q = 0.1,
  hd_thresh = 0)

umap_two_curvy_scaled <- algo_obj_two_curvy$nldr_scaled_obj$scaled_nldr
tr_from_to_df_two_curvy <- algo_obj_two_curvy$trimesh_data
df_bin_centroids_two_curvy <- algo_obj_two_curvy$model_2d
df_bin_two_curvy <- algo_obj_two_curvy$model_highd
hex_grid <- algo_obj_two_curvy$hb_obj$hex_poly
counts_df <- algo_obj_two_curvy$hb_obj$std_cts

## Predict

predict_df3 <- predict_emb(highd_data = test_data_two_curvy[,1:8],
                           model_highd = df_bin_two_curvy,
                           model_2d = df_bin_centroids_two_curvy)

## Approach A: Compare train vs test UMAP positions

ggplot(training_umap_two_curvy,
       aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha = 0.5, colour = "#d0d1e6", size = 0.5) +
  geom_point(data = test_umap_two_curvy, aes(
    x = UMAP1, y = UMAP2
  ), alpha = 0.5, colour = "#045a8d", size = 0.5)


## Approach B: Compare UMAP true test positions vs predicted positions

ggplot(test_umap_two_curvy,
       aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha = 0.5, colour = "#d0d1e6", size = 0.5) +
  geom_point(data = predict_UMAP_df, aes(
    x = UMAP1, y = UMAP2
  ), alpha = 0.5, colour = "#045a8d", size = 0.5)


## Approach C: Compare quollr positions vs UMAP positions

quollr_predict_umap_test_scaled <- gen_scaled_data(nldr_data = test_umap_two_curvy)$scaled_nldr

ggplot(quollr_predict_umap_test_scaled,
       aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha = 0.5, colour = "#d0d1e6", size = 0.5) +
  geom_point(data = predict_df3, aes(
    x = pred_emb_1, y = pred_emb_2
  ), alpha = 0.5, colour = "#045a8d", size = 0.5)


################################################################################

# tSNE simulate new test data and then compare training vs new predicted test



# PHATE simulate new test data and then compare training vs new predicted test


################################################################################

# New badly simulate outside domain of training and compare for best tSNE
