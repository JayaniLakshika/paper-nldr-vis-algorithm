library(tibble)
library(dplyr)
library(langevitour)
library(readr)

set.seed(20240110)

data <- read_rds(here::here("data/two_c_shaped_clust/two_c_shaped_clust_data.rds"))
model_data <- read_rds(here::here("data/two_c_shaped_clust/two_c_shaped_clust_true_model.rds"))
connections <- read_rds(here::here("data/two_c_shaped_clust/two_c_shaped_clust_true_model_connections.rds"))

df <- bind_rows(model_data |> select(-ID) |> mutate(type = "model"),
                data |> mutate(type = "data"))

langevitour(df |> select(-type),
            group = df$type,
            lineFrom = connections$from,
            lineTo = connections$to)
