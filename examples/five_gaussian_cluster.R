library(dplyr)
library(snedata)
library(ggflowchart)
library(purrr) ## map function
library(gridExtra) ## for grid.arrange
library(rsample)
library(DT)
library(ggbeeswarm)
library(ggplot2)
library(readr)

library(Rtsne)
library(umap)
library(phateR)
library(reticulate)
library(patchwork)

library(grid)

set.seed(20240110)

source("quollr_code.R", local = TRUE)
source("nldr_code.R", local = TRUE)

sample_size <- 5000
cluster_size <- sample_size/5
df1 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.05), y=rnorm(cluster_size, mean = 0, sd = 0.05), z=rnorm(cluster_size, mean = 0, sd = 0.05), w=rnorm(cluster_size, mean = 0, sd = 0.05))

df2 <- tibble::tibble(x=rnorm(cluster_size, mean = 1, sd = 0.05), y=rnorm(cluster_size, mean = 0, sd = 0.05), z=rnorm(cluster_size, mean = 0, sd = 0.05), w=rnorm(cluster_size, mean = 0, sd = 0.05))

df3 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.05), y=rnorm(cluster_size, mean = 1, sd = 0.05), z=rnorm(cluster_size, mean = 0, sd = 0.05), w=rnorm(cluster_size, mean = 0, sd = 0.05))

df4 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.05), y=rnorm(cluster_size, mean = 0, sd = 0.05), z=rnorm(cluster_size, mean = 1, sd = 0.05), w=rnorm(cluster_size, mean = 0, sd = 0.05))

df5 <- tibble::tibble(x=rnorm(cluster_size, mean = 0, sd = 0.05), y=rnorm(cluster_size, mean = 0, sd = 0.05), z=rnorm(cluster_size, mean = 0, sd = 0.05), w=rnorm(cluster_size, mean = 1, sd = 0.05))

df_2 <- bind_rows(df1, df2, df3, df4, df5)
df_2 <- df_2 |>
  rename(x1 = x, x2 = y, x3 = z, x4 = w)

df_2 <- df_2 |>
  mutate(ID = row_number())

data_split_sp <- initial_split(df_2)
training_data_5 <- training(data_split_sp) |>
  arrange(ID)
test_data_5 <- testing(data_split_sp) |>
  arrange(ID)

### tSNE
tSNE_data_gau <- Fit_tSNE(training_data_5 |> dplyr::select(-ID), opt_perplexity = calculate_effective_perplexity(training_data_5), with_seed = 20240110)

plot_gau <- plot_tSNE_2D(tSNE_data_gau) + #ggtitle("(a)") +
  theme_linedraw() +
  theme(plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### UMAP

UMAP_fit <- umap(training_data_5 |> dplyr::select(-ID), n_neighbors = 15, n_components =  2)

UMAP_data_gau <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data_gau)[1:(ncol(UMAP_data_gau))] <- paste0(rep("UMAP",(ncol(UMAP_data_gau))), 1:(ncol(UMAP_data_gau)))

UMAP_data_gau <- UMAP_data_gau |>
  mutate(ID = training_data_5$ID)

## predict umap embeddings

predict_UMAP_df <- predict(UMAP_fit, test_data_5 |> dplyr::select(-ID)) |>
  as.data.frame()

names(predict_UMAP_df)[1:(ncol(predict_UMAP_df))] <- paste0(rep("UMAP",(ncol(predict_UMAP_df))), 1:(ncol(predict_UMAP_df)))

predict_UMAP_df <- predict_UMAP_df |>
  mutate(ID = test_data_5$ID)

plot_UMAP_2D(UMAP_data_gau) +
  geom_point(data = predict_UMAP_df, aes(x = UMAP1, y = UMAP2), color = "red")


### TriMAP


### PaCMAP


### Phate

PHATE_data <- Fit_PHATE(training_data_5, knn = 5, with_seed = 20240110)
write_csv(PHATE_data, paste0(here::here(), "/data/phate_data_s_curve.csv"))
