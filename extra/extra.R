```{r}
#| echo: false
#| warning: false

error_df_lg_1 <- read_rds("data/s_curve/s_curve_summary_lg_tsne.rds")
error_df_lg_2 <- read_rds("data/s_curve/s_curve_summary_lg_umap.rds")
error_df_lg_3 <- read_rds("data/s_curve/s_curve_summary_lg_phate.rds")
error_df_lg_4 <- read_rds("data/s_curve/s_curve_summary_lg_trimap.rds")
error_df_lg_5 <- read_rds("data/s_curve/s_curve_summary_lg_pacmap.rds")

error_df_lg <- dplyr::bind_rows(error_df_lg_1, error_df_lg_2, error_df_lg_3, error_df_lg_4, error_df_lg_5)

error_df_lg$method <- factor(error_df_lg$method, levels = c("tSNE", "UMAP", "PHATE", "TriMAP", "PaCMAP"))

## To draw with AIC
aic_plot_lg <- ggplot(error_df_lg |> dplyr::filter(method == "UMAP"), aes(x = benchmark_rm_lg,
                                                                          y = total_error,
                                                                          color = method
)) +
  geom_point() +
  geom_vline(xintercept = 1.347, linetype="dashed",
             color = "red", size=0.5) +
  geom_line() +
  #geom_vline(xintercept = NROW(full_grid_with_hexbin_id)) +
  #annotate("text", x= (NROW(full_grid_with_hexbin_id) - 10), y=-5000, label=paste0("effective number of bins = ", as.character(NROW(full_grid_with_hexbin_id))), angle=90) +
  theme_light() +
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7)) +
  scale_colour_manual(values = c("#377eb8", "#e41a1c", "#4daf4a", "#984ea3", "#ff7f00")) +
  ylab("AIC") +
  xlab("benchmark value")
## Effective number of bins along x-axis

mse_plot_lg <- ggplot(error_df_lg |> dplyr::filter(method == "UMAP"), aes(x = benchmark_rm_lg,
                                                                          y = total_mse,
                                                                          color = method
)) +
  geom_point() +
  geom_vline(xintercept = 1.347, linetype="dashed",
             color = "red", size=0.5) +
  geom_line() +
  theme_light() +
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7)) +
  # geom_vline(xintercept = NROW(full_grid_with_hexbin_id)) +
  # annotate("text", x= (NROW(full_grid_with_hexbin_id) - 10), y=0.25, label=paste0("effective number of bins = ", as.character(NROW(full_grid_with_hexbin_id))), angle=90) +
  scale_colour_manual(values = c("#377eb8", "#e41a1c", "#4daf4a", "#984ea3", "#ff7f00")) +
  ylab("MSE") +
  xlab("benchmark value")

```

```{r}
#| echo: false
#| warning: false
#| fig-cap: "Goodness of fit statistics from UMAP applied to training S-curve dataset with different benchmark values to remove the long edges. What is the effective benchmark value to remove the long edges? $1.347$ is the chosen benchmark value to remove long edges. The value is selected through a visual inspection of the triangular mesh is conducted to verify whether it successfully eliminates long edges."
#| label: fig-diagnosticpltScurvelgrm
##| out-width: 100%
#| fig-pos: H

aic_plot_lg + mse_plot_lg +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides='collect', ncol = 2) &
  theme(legend.position='none', plot.tag = element_text(size = 8))
```
