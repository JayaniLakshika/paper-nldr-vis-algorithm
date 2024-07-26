library(crosstalk)
library(plotly)
library(readr)
library(quollr)
library(langevitour)

umap_scurve <- read_rds(file = "data/s_curve/s_curve_umap.rds")

scurve_scaled_obj <- gen_scaled_data(
  data = umap_scurve)

umap_scurve_scaled <- scurve_scaled_obj$scaled_nldr
lim1 <- scurve_scaled_obj$lim1
lim2 <- scurve_scaled_obj$lim2
r2 <- diff(lim2)/diff(lim1)

df_all_scurve <- dplyr::bind_cols(training_data_scurve |> dplyr::select(-ID),
                                  umap_data_with_hb_id)

# df_b <- df_bin_scurve |>
#   dplyr::filter(hb_id %in% df_bin_centroids_scurve$hexID) |>
#   dplyr::mutate(type = "model") ## Data with summarized mean
#
# ## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
# df_b_n <- df_b[match(df_bin_centroids_scurve$hexID, df_b$hb_id),] |>
#   select(-type)
#
# df_b_n <- left_join(df_b_n, df_bin_centroids_scurve, by = c("hb_id" = "hexID")) |>
#   select(-hb_id, -std_counts, -drop_empty)

df_all_scurve_n <- df_all_scurve |>
  select(-ID, -hb_id)

# names(df_b_n) <- names(df_all_scurve_n)

shared_df_scurve <- SharedData$new(df_all_scurve_n)

nldr_scurve <- shared_df_scurve |>
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha=0.5, colour="#6a3d9a", size = 0.5) +
  theme_linedraw() +
  theme(
    #aspect.ratio = 1,
    plot.background = element_rect(fill = 'transparent', colour = NA),
    plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
    panel.background = element_rect(fill = 'transparent',
                                    colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )

nldr_scurve_plt <- ggplotly(nldr_scurve, width = as.character(round(600/r2, 0)),
                            height = "600") |>
  style(unselected=list(marker=list(opacity=1))) |>
  highlight(on="plotly_selected", off="plotly_deselect") |>
  config(displayModeBar = FALSE)

langevitour_output <- langevitour(training_data |> select(-type),
                                  link=shared_df_scurve)

bscols(nldr_scurve_plt, langevitour_output)
