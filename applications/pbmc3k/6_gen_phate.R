library(phateR)
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

PHATE_data <- Fit_PHATE(training_data_pbmc |> dplyr::select(-ID), knn = 10, with_seed = 20240110)
PHATE_data <- PHATE_data |>
  select(PHATE1, PHATE2)
PHATE_data <- PHATE_data |>
  mutate(ID = training_data_pbmc$ID)

write_rds(PHATE_data, "data/pbmc3k/pbmc_phate_10.rds")

phate_pbmc <- PHATE_data
phate_pbmc_scaled <- as.data.frame(do.call(cbind, gen_scaled_data(data = phate_pbmc, x = "PHATE1", y = "PHATE2"))) |>
  dplyr::rename(c("PHATE1" = "scaled_PHATE1",
                  "PHATE2" = "scaled_PHATE2")) |>
  dplyr::mutate(ID = 1:NROW(phate_pbmc))

phate_plot_pbmc <- phate_pbmc_scaled |>
  ggplot(aes(x = PHATE1,
             y = PHATE2))+
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

phate_plot_pbmc
