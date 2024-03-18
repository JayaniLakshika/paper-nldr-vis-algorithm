## Import data
training_data_pbmc <- read_rds("data/pbmc/pbmc_3k_festem/pbmc_pca_50.rds")
training_data_pbmc <- training_data_pbmc[, 1:15] |>
  mutate(ID = 1:NROW(training_data_pbmc))

umap_pbmc <- read_rds("data/pbmc/pbmc_umap_5_min_dist_0.99_metric_cosine.rds") |>
  dplyr::select(-c(ID, cell_label))

dist_vec <- proxy::dist(x = umap_pbmc, method = "Euclidean") |> as.vector()

from_vec <- c()
to_vec <- c()
num_obs <- 1:(NROW(umap_pbmc) - 1)

for (obs in num_obs) {

  from_val <- rep(obs, (NROW(umap_pbmc) - obs))
  if ((obs + 1) <= NROW(umap_pbmc)) {
    to_val <- (obs + 1):NROW(umap_pbmc)
  }
  from_vec <- append(from_vec, from_val)
  to_vec <- append(to_vec, to_val)

}

dist_2d <- tibble::tibble(from = from_vec, to = to_vec, dist = dist_vec)

## Import detour results
detour_results8 <- read_csv("data/pbmc/pbmc_3k_festem/detourr_export8.csv")

### Rename labels

detour_results8 <- detour_results8 |>
  dplyr::mutate(result2_class_label = case_when(
    colour == "000000" ~ "cluster 4",
    colour == "49df20" ~ "cluster 2",
    colour == "619cff" ~ "cluster 3",
    colour == "d82cd3" ~ "cluster 1",
    colour == "df2020" ~ "cluster 5",
    colour == "ffbd61" ~ "cluster 6",
    .default = "no cluster"

  ))

detour_results8_n <- detour_results8 |>
  dplyr::mutate(ID = training_data_pbmc$ID) |>
  dplyr::select(ID, result2_class_label)

## Map cluster labels for from and to
dist_2d <- dplyr::inner_join(dist_2d, detour_results8_n, by = c("from" = "ID")) |>
  dplyr::rename("from_cluster" = "result2_class_label")

dist_2d <- dplyr::inner_join(dist_2d, detour_results8_n, by = c("to" = "ID")) |>
  dplyr::rename("to_cluster" = "result2_class_label")

write_rds(dist_2d, "data/pbmc/pbmc_3k_festem/lowd_dist_pbmc.rds")

### Compute high-D distances

dist_vec <- proxy::dist(x = training_data_pbmc[, -16], method = "Euclidean") |> as.vector()

from_vec <- c()
to_vec <- c()
num_obs <- 1:(NROW(training_data_pbmc) - 1)

for (obs in num_obs) {

  from_val <- rep(obs, (NROW(training_data_pbmc) - obs))
  if ((obs + 1) <= NROW(training_data_pbmc)) {
    to_val <- (obs + 1):NROW(training_data_pbmc)
  }
  from_vec <- append(from_vec, from_val)
  to_vec <- append(to_vec, to_val)

}

dist_highd <- tibble::tibble(from = from_vec, to = to_vec, dist = dist_vec)


## Map cluster labels for from and to
dist_highd <- dplyr::inner_join(dist_highd, detour_results8_n, by = c("from" = "ID")) |>
  dplyr::rename("from_cluster" = "result2_class_label")

dist_highd <- dplyr::inner_join(dist_highd, detour_results8_n, by = c("to" = "ID")) |>
  dplyr::rename("to_cluster" = "result2_class_label")

write_rds(dist_highd, "data/pbmc/pbmc_3k_festem/highd_dist_pbmc.rds")


umap_pbmc <- read_rds("data/pbmc/pbmc_3k_festem/pbmc_umap.rds") |>
  dplyr::select(-cell_label)

dist_vec <- proxy::dist(x = umap_pbmc, method = "Euclidean") |> as.vector()

from_vec <- c()
to_vec <- c()
num_obs <- 1:(NROW(umap_pbmc) - 1)

for (obs in num_obs) {

  from_val <- rep(obs, (NROW(umap_pbmc) - obs))
  if ((obs + 1) <= NROW(umap_pbmc)) {
    to_val <- (obs + 1):NROW(umap_pbmc)
  }
  from_vec <- append(from_vec, from_val)
  to_vec <- append(to_vec, to_val)

}

dist_2d <- tibble::tibble(from = from_vec, to = to_vec, dist = dist_vec)
write_rds(dist_2d, "data/pbmc/pbmc_3k_festem/lowd_dist_pbmc_param1.rds")

dist_2d <- dist_2d |>
  dplyr::rename("dist_2d" = "dist")

dist_highd <- dist_highd |>
  dplyr::rename("dist_highd" = "dist") |>
  dplyr::select(dist_highd)

dist_df <- dplyr::bind_cols(dist_2d, dist_highd)
write_rds(dist_df, "data/mnist/dist_pbmc_param1.rds")
