library(tidyverse)
data <- read_rds(here::here("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_data.rds"))
tSNE_data <- read_rds(here::here("data/two_nonlinear/two_non_linear_diff_shaped_close_clusters_tsne_perplexity_47.rds"))


score_result <- my_instance$global_score(as.matrix(data), as.matrix(tSNE_data))
score_result
results <- py$evaluate_output(X = as.matrix(data), X_new = as.matrix(tSNE_data), y = sample(0:4, 2000, replace = TRUE), name = "my_run", baseline = FALSE, labelled = TRUE)
results
