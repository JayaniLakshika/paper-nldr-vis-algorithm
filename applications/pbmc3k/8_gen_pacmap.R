library(reticulate)
library(readr)
library(dplyr)

set.seed(20240110)

use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

reticulate::source_python(paste0(here::here(), "/examples/function_scripts/Fit_PacMAP_code.py"))
reticulate::source_python(paste0(here::here(), "/examples/function_scripts/Fit_TriMAP_code.py"))

source("nldr_code.R", local = TRUE)

## Import data
training_data_pbmc <- read_rds("data/pbmc3k/pbmc_pca_50.rds")
training_data_pbmc <- training_data_pbmc[, 1:9] |>
  mutate(ID = 1:NROW(training_data_pbmc))

## PaCMAP

tem_dir <- tempdir()

Fit_PacMAP_data(training_data_pbmc |> dplyr::select(-ID), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_PaCMAP_values.csv")

Fit_PaCMAP(as.integer(2), as.integer(30), "random", 0.1, as.integer(1), path, path2)

PacMAP_data <- read_csv(path2)
PacMAP_data <- PacMAP_data |>
  mutate(ID = training_data_pbmc$ID)

write_rds(PacMAP_data, file = "data/pbmc3k/pbmc_pacmap_30_random_0.1_1.rds")

pacmap_pbmc <- PacMAP_data
pacmap_pbmc_scaled <- as.data.frame(do.call(cbind, gen_scaled_data(data = pacmap_pbmc, x = "PaCMAP1", y = "PaCMAP2"))) |>
  dplyr::rename(c("PaCMAP1" = "scaled_PaCMAP1",
                  "PaCMAP2" = "scaled_PaCMAP2")) |>
  dplyr::mutate(ID = 1:NROW(pacmap_pbmc))

pacmap_plot_pbmc <- pacmap_pbmc_scaled |>
  ggplot(aes(x = PaCMAP1,
             y = PaCMAP2))+
  geom_point(alpha=0.5) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
  theme_linedraw() +
  theme(legend.position = "none", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'))

pacmap_plot_pbmc
