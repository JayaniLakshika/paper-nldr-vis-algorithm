library(readr)

### Spin-brush
#remotes::install_github("casperhart/detourr")
library(detourr)

detour(training_data_pbmc,
       tour_aes(projection = PC_1:PC_15)) |>
  tour_path(grand_tour(2), fps = 60,
            max_bases=20) |>
  show_scatter(alpha = 0.7,
               axes = FALSE)

## Import detour results
detour_results1 <- read_csv("data/pbmc/pbmc_3k_festem/detourr_export.csv")
#names(detour_results)[1:15] <- paste0("PC_", 1:15)
table(training_data_pbmc$cell_label, detour_results1$colour)

## Import detour results
detour_results2 <- read_csv("data/pbmc/pbmc_3k_festem/detourr_export2.csv")
#names(detour_results)[1:15] <- paste0("PC_", 1:15)
table(training_data_pbmc$cell_label, detour_results2$colour)

## Import detour results
detour_results3 <- read_csv("data/pbmc/pbmc_3k_festem/detourr_export3.csv")
#names(detour_results)[1:15] <- paste0("PC_", 1:15)
table(training_data_pbmc$cell_label, detour_results3$colour)


### Spin-brush with detour and model
detour(df_exe,
       tour_aes(projection = PC_1:PC_15,
                colour = type)) |>
  tour_path(grand_tour(2), fps = 60,
            max_bases=20) |>
  show_scatter(alpha = 0.7,
               axes = FALSE,
               edges = distance_df_small_edges |> select(from, to) |> as.matrix(),
               palette = c("#fee0d2", "#fcbba1",
                           "#fc9272", "#fb6a4a", "#ef3b2c",
                           "#cb181d", "#a50f15", "#99000d", "#33a02c"))


### Experiment with liminal
library(liminal)

training_data_pbmc <- training_data_pbmc |>
  mutate(cell_label = UMAP_pbmc$cell_label)

limn_tour_link(embed_data = UMAP_pbmc |> select("UMAP1", "UMAP2"),
               tour_data = training_data_pbmc,
               cols = PC_1:PC_15, # tour columns to select
               color = cell_label # variable to highlight across both view, can come for either data frames
)
