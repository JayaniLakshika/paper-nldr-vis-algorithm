library(readr)
library(quollr)
library(dplyr)

data_mnist <- read_rds("data/mnist/mnist_10_pcs_of_digit_1.rds")
names(data_mnist) <- paste0("x", 1:NCOL(data_mnist))

data_mnist <- data_mnist |>
  mutate(ID = 1:NROW(data_mnist))

## For umap
umap_mnist <- read_rds("data/mnist/mnist_umap.rds")
error_mnist_umap <- gen_diffbin1_errors(highd_data = data_mnist, nldr_data = umap_mnist) |>
  dplyr::mutate(method = "UMAP")

write_rds(error_mnist_umap, "data/mnist/error_mnist_umap.rds")

###########

## For tsne
#tsne_mnist <- read_rds("data/mnist/mnist_tsne89.rds")
tsne_mnist <- read_rds("data/mnist/mnist_tsne30.rds")

error_mnist_tsne <- gen_diffbin1_errors(highd_data = data_mnist, nldr_data = tsne_mnist) |>
  dplyr::mutate(method = "tSNE")

write_rds(error_mnist_tsne, "data/mnist/error_mnist_tsne.rds")

###########

## For tsne
tsne_mnist <- read_rds("data/mnist/mnist_tsne89.rds")
error_mnist_tsne <- gen_diffbin1_errors(highd_data = data_mnist, nldr_data = tsne_mnist) |>
  dplyr::mutate(method = "tSNE2")

write_rds(error_mnist_tsne, "data/mnist/error_mnist_tsne2.rds")


###########


## For phate
phate_mnist <- read_rds("data/mnist/mnist_phate.rds")
error_mnist_phate <- gen_diffbin1_errors(highd_data = data_mnist, nldr_data = phate_mnist) |>
  dplyr::mutate(method = "PHATE")

write_rds(error_mnist_phate, "data/mnist/error_mnist_phate.rds")

###########

## For trimap
trimap_mnist <- read_rds("data/mnist/mnist_trimap.rds")
error_mnist_trimap <- gen_diffbin1_errors(highd_data = data_mnist, nldr_data = trimap_mnist) |>
  dplyr::mutate(method = "TriMAP")

write_rds(error_mnist_trimap, "data/mnist/error_mnist_trimap.rds")

###########
## For pacmap
pacmap_mnist <- read_rds("data/mnist/mnist_pacmap.rds")
error_mnist_pacmap <- gen_diffbin1_errors(highd_data = data_mnist, nldr_data = pacmap_mnist) |>
  dplyr::mutate(method = "PaCMAP")

write_rds(error_mnist_pacmap, "data/mnist/error_mnist_pacmap.rds")
