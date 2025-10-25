## This script is to diagnose layout e of MNIST data

library(quollr)
library(tidyverse)
set.seed(20240110)

## Data
data_mnist <- read_rds("data/mnist/mnist_10_pcs_of_digit_1.rds")
names(data_mnist) <- paste0("x", 1:NCOL(data_mnist))

data_mnist <- data_mnist |>
  mutate(ID = 1:NROW(data_mnist))

pacmap_mnist <- read_rds("data/mnist/mnist_pacmap.rds")

mnist_pacmap_obj <- fit_highd_model(highd_data = data_mnist,
                                    nldr_data = pacmap_mnist,
                                    b1 = 30,
                                    q = 0.1,
                                    hd_thresh = 0)

model_error <- augment(x = mnist_pacmap_obj,
                       highd_data = data_mnist)

df_exe <- comb_all_data_model_error(highd_data = data_mnist,
                                    nldr_data = pacmap_mnist,
                                    model_highd = mnist_pacmap_obj$model_highd,
                                    model_2d = mnist_pacmap_obj$model_2d,
                                    error_data = model_error)

edge_data <- scurve_model_obj$trimesh_data

show_error_link_plots(point_data = df_exe, edge_data = edge_data)
