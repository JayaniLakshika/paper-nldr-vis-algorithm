library(readr)
library(quollr)
library(dplyr)

set.seed(20240110)

quad <- function(a = 3, b = 2 * a2, c = -(a2^2 + a1^2))
{
  a <- as.complex(a)
  answer <- c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
              (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
  if(all(Im(answer) == 0)) answer <- Re(answer)
  if(answer[1] == answer[2]) return(answer[1])
  answer[answer>0] ## only positive
}


data_two_curvy <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_data.rds")
data_two_curvy <- data_two_curvy |>
  mutate(ID = 1:NROW(data_two_curvy))

## For umap
umap_two_non_linear_diff_shaped_close_clusters <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_umap_n-neigbors_15_min-dist_0.1.rds") |>
  mutate(ID = row_number())

error_two_non_linear_diff_shaped_close_clusters_umap <- gen_diffbin1_errors(highd_data = data_two_curvy,
                                                                            nldr_data = umap_two_non_linear_diff_shaped_close_clusters) |>
  dplyr::mutate(side_length = quad(a=3, b = 2 * a2, c = -(a2^2 + a1^2))) |>
  dplyr::mutate(method = "UMAP")

write_rds(error_two_non_linear_diff_shaped_close_clusters_umap, "data/two_nonlinear/error_two_non_linear_diff_shaped_close_clusters_umap.rds")

###########

## For tsne
#tsne_two_non_linear_diff_shaped_close_clusters <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_tsne89.rds")
tsne_two_non_linear_diff_shaped_close_clusters <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_47.rds") |>
  mutate(ID = row_number())

error_two_non_linear_diff_shaped_close_clusters_tsne <- gen_diffbin1_errors(highd_data = data_two_curvy,
                                                                            nldr_data = tsne_two_non_linear_diff_shaped_close_clusters) |>
  dplyr::mutate(side_length = quad(a=3, b = 2 * a2, c = -(a2^2 + a1^2))) |>
  dplyr::mutate(method = "tSNE")

write_rds(error_two_non_linear_diff_shaped_close_clusters_tsne, "data/two_nonlinear/error_two_non_linear_diff_shaped_close_clusters_tsne.rds")

###########

## For tsne
tsne_two_non_linear_diff_shaped_close_clusters <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_62.rds") |>
  mutate(ID = row_number())

error_two_non_linear_diff_shaped_close_clusters_tsne <- gen_diffbin1_errors(highd_data = data_two_curvy,
                                                                            nldr_data = tsne_two_non_linear_diff_shaped_close_clusters) |>
  dplyr::mutate(side_length = quad(a=3, b = 2 * a2, c = -(a2^2 + a1^2))) |>
  dplyr::mutate(method = "tSNE2")

write_rds(error_two_non_linear_diff_shaped_close_clusters_tsne, "data/two_nonlinear/error_two_non_linear_diff_shaped_close_clusters_tsne2.rds")


###########

## For phate
phate_two_non_linear_diff_shaped_close_clusters <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_phate_knn_5.rds") |>
  mutate(ID = row_number())

error_two_non_linear_diff_shaped_close_clusters_phate <- gen_diffbin1_errors(highd_data = data_two_curvy,
                                                                            nldr_data = phate_two_non_linear_diff_shaped_close_clusters) |>
  dplyr::mutate(side_length = quad(a=3, b = 2 * a2, c = -(a2^2 + a1^2))) |>
  dplyr::mutate(method = "PHATE")

write_rds(error_two_non_linear_diff_shaped_close_clusters_phate, "data/two_nonlinear/error_two_non_linear_diff_shaped_close_clusters_phate.rds")

###########

## For trimap
trimap_two_non_linear_diff_shaped_close_clusters <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_trimap_n-inliers_12_n-outliers_4_n-random_3.rds") |>
  mutate(ID = row_number())

error_two_non_linear_diff_shaped_close_clusters_trimap <- gen_diffbin1_errors(highd_data = data_two_curvy,
                                                                             nldr_data = trimap_two_non_linear_diff_shaped_close_clusters) |>
  dplyr::mutate(side_length = quad(a=3, b = 2 * a2, c = -(a2^2 + a1^2))) |>
  dplyr::mutate(method = "TriMAP")

write_rds(error_two_non_linear_diff_shaped_close_clusters_trimap, "data/two_nonlinear/error_two_non_linear_diff_shaped_close_clusters_trimap.rds")

###########

## For pacmap
pacmap_two_non_linear_diff_shaped_close_clusters <- read_rds("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_pacmap_n-neighbors_10_init_random_MN-ratio_0.5_FP-ratio_2.rds") |>
  mutate(ID = row_number())

error_two_non_linear_diff_shaped_close_clusters_pacmap <- gen_diffbin1_errors(highd_data = data_two_curvy,
                                                                              nldr_data = pacmap_two_non_linear_diff_shaped_close_clusters) |>
  dplyr::mutate(side_length = quad(a=3, b = 2 * a2, c = -(a2^2 + a1^2))) |>
  dplyr::mutate(method = "PaCMAP")

write_rds(error_two_non_linear_diff_shaped_close_clusters_pacmap, "data/two_nonlinear/error_two_non_linear_diff_shaped_close_clusters_pacmap.rds")
