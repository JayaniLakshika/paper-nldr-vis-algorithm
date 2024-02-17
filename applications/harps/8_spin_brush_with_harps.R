library(readr)
library(dplyr)

### Spin-brush
#remotes::install_github("casperhart/detourr")
library(detourr)

training_data_harps_with_label <- read_rds("data/harps/harps_data.rds")
#UMAP_harps <- read_rds("data/pbmc/pbmc_3k_festem/pbmc_umap.rds")

training_data_harps <- training_data_harps_with_label[, 1:17] |>
  mutate(ID = 1:NROW(training_data_harps_with_label))

names(training_data_harps) <- append(paste0("x", 1:17), "ID")

detour(training_data_harps,
       tour_aes(projection = x1:x17)) |>
  tour_path(grand_tour(2), fps = 30,
            max_bases=100) |>
  show_scatter(alpha = 0.7,
               axes = FALSE)


## Import detour results
detour_results1 <- read_csv("data/harps/detourr_export1.csv")
detour_results1 <- detour_results1 |>
  dplyr::mutate(ID = training_data_harps$ID)

table(training_data_harps_with_label$cluster, detour_results1$colour)
