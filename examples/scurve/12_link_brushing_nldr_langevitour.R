df_all_scurve_n <- df_all_scurve |>
  select(-ID, -hb_id) |>
  dplyr::mutate(type = "data")


df_b <- df_bin_scurve |>
  dplyr::filter(hb_id %in% df_bin_centroids_scurve$hexID)

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_scurve$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id) |>
  dplyr::mutate(type = "model")

# df_bin_centroids_scurve_n <-df_bin_centroids_scurve[match(df_b$hb_id, df_bin_centroids_scurve$hexID),] |>
#   select(c_x, c_y) |>
#   rename(c("UMAP1" = "c_x",
#            "UMAP2" = "c_y"))
#
# df_b <- bind_cols(df_b, df_bin_centroids_scurve_n) |>
#   dplyr::select(-hb_id) |>
#   dplyr::mutate(type = "model")

df_exe <- dplyr::bind_rows(df_b, df_all_scurve_n)

## Set the maximum difference as the criteria
distance_df_small_edges <- distance_scurve |>
  dplyr::filter(distance < benchmark_scurve)

shared_df_scurve <- SharedData$new(df_exe)


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


langevitour_output <- langevitour::langevitour(df_exe[1:7],
                                               lineFrom = distance_df_small_edges$from,
                                               lineTo = distance_df_small_edges$to,
                                               group = df_exe$type, pointSize = append(rep(2, NROW(df_b)), rep(1, NROW(df))),
                                               levelColors = c("#6a3d9a", "#33a02c"),
                                               link=shared_df_scurve,
                                               link_filter=TRUE)

bscols(nldr_scurve_plt, langevitour_output)
