## Import data
training_data_mnist <- read_rds("data/mnist/mnist_10_pcs_of_digit_1.rds")
training_data_mnist <- training_data_mnist |>
  mutate(ID = 1:NROW(training_data_mnist))

pacmap_minst <- read_rds("data/mnist/mnist_pacmap.rds")

dist_vec <- proxy::dist(x = pacmap_minst, method = "Euclidean") |> as.vector()

from_vec <- c()
to_vec <- c()
num_obs <- 1:(NROW(pacmap_minst) - 1)

for (obs in num_obs) {

  from_val <- rep(obs, (NROW(pacmap_minst) - obs))
  if ((obs + 1) <= NROW(pacmap_minst)) {
    to_val <- (obs + 1):NROW(pacmap_minst)
  }
  from_vec <- append(from_vec, from_val)
  to_vec <- append(to_vec, to_val)

}

dist_2d <- tibble::tibble(from = from_vec, to = to_vec, dist = dist_vec)

dist_2d <- dist_2d |>
  dplyr::rename("dist_2d" = "dist")

write_rds(dist_2d, "data/mnist/lowd_dist_mnist.rds")

### Compute high-D distances

dist_vec <- proxy::dist(x = training_data_mnist[, -11], method = "Euclidean") |> as.vector()

from_vec <- c()
to_vec <- c()
num_obs <- 1:(NROW(training_data_mnist) - 1)

for (obs in num_obs) {

  from_val <- rep(obs, (NROW(training_data_mnist) - obs))
  if ((obs + 1) <= NROW(training_data_mnist)) {
    to_val <- (obs + 1):NROW(training_data_mnist)
  }
  from_vec <- append(from_vec, from_val)
  to_vec <- append(to_vec, to_val)

}

dist_highd <- tibble::tibble(from = from_vec, to = to_vec, dist = dist_vec)

dist_highd <- dist_highd |>
  dplyr::rename("dist_highd" = "dist") |>
  dplyr::select(dist_highd)

write_rds(dist_highd, "data/mnist/highd_dist_mnist.rds")

dist_df <- dplyr::bind_cols(dist_2d, dist_highd)
write_rds(dist_df, "data/mnist/dist_mnist.rds")
