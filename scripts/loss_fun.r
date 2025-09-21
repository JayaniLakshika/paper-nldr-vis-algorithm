library(dplyr)
library(tibble)
library(readr)
library(conflicted)

library(Rtsne)
library(uwot)

set.seed(20240110)

data <- read_rds(here::here("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_data.rds"))

## tSNE (default)
perplexity <- 47

tSNE_fit <- data |>
  dplyr::select(where(is.numeric)) |>
  Rtsne::Rtsne(perplexity = perplexity,
               pca = FALSE)

## Loss function value
tail(tSNE_fit$itercosts, 1)

## UMAP (default)
### https://github.com/sciai-lab/UMAPs-true-loss

losses <- numeric()

callback <- function(epoch, n_epochs, coords) {
  # Simple placeholder: you could compute a UMAP-like loss here
  # For demonstration, let's store the mean squared distance from origin
  losses[epoch] <<- mean(rowSums(coords^2))
  cat(sprintf("Epoch %d/%d finished. Loss approx: %.4f\n", epoch, n_epochs, losses[epoch]))
}

umap_fit <- uwot::umap(
  X = data,
  n_neighbors = 15,
  n_components = 2,
  n_epochs = 100,
  verbose = TRUE,
  ret_model = TRUE,
  epoch_callback = callback
)

## Loss function value
tail(losses, 1)






