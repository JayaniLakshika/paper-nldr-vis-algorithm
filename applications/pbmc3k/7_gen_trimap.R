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

## TriMAP

tem_dir <- tempdir()

Fit_TriMAP_data(training_data_pbmc |> dplyr::select(-ID), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_TriMAP_values.csv")

Fit_TriMAP(as.integer(2), as.integer(30), as.integer(4), as.integer(10), path, path2)

TriMAP_data <- read_csv(path2)
TriMAP_data <- TriMAP_data |>
  mutate(ID = training_data_pbmc$ID)

write_rds(TriMAP_data, file = "data/pbmc3k/pbmc_trimap_30_4_10.rds")

trimap_pbmc <- TriMAP_data
trimap_pbmc_scaled <- as.data.frame(do.call(cbind, gen_scaled_data(data = trimap_pbmc, x = "TriMAP1", y = "TriMAP2"))) |>
  dplyr::rename(c("TriMAP1" = "scaled_TriMAP1",
                  "TriMAP2" = "scaled_TriMAP2")) |>
  dplyr::mutate(ID = 1:NROW(trimap_pbmc))

trimap_plot_pbmc <- trimap_pbmc_scaled |>
  ggplot(aes(x = TriMAP1,
             y = TriMAP2))+
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

trimap_plot_pbmc
