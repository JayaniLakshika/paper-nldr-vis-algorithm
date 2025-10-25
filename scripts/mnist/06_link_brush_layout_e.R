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

show_error_link_plots(point_data = df_exe, edge_data = edge_data, point_colours = c("red", "#FF7755"))

## Let's dig more
## Image ID's related to the outlier points: 3619, 3554, 4678, 3899, 3486, 5975, 501, 1012, 2260, 3352, 5747, 27779

library(plotly)
ggplot(pacmap_mnist, aes(x = emb1, y = emb2, label = ID)) +
  geom_text() +
  theme(aspect.ratio = 1)

ggplotly()

source("scripts/additional_functions.R")

## Data with pixel values
mnist_data <- read_rds("data/mnist/mnist_digit_1.rds")

mnist_data <- mnist_data |>
  mutate(instance = row_number()) |>
  gather(pixel, value, -Label, -instance) |>
  extract(pixel, "pixel", "(\\d+)", convert = TRUE) |>
  mutate(pixel = pixel - 2, x = pixel %% 28, y = 28 - pixel %/% 28)

img_outliers <- c(3619, 3554, 4678, 3899, 3486, 5975, 501, 1012, 2260, 3352, 5747, 27779)

pixels_gathered_out <-  mnist_data |>
  filter(instance %in% img_outliers)

plot_digit_img(
  pixels_gathered_out, color_on = "#375cb1", title_text = "b")
