library(crosstalk)
library(plotly)

umap_scurve <- read_rds(file = "data/s_curve/s_curve_umap.rds")

scurve_scaled_obj <- gen_scaled_data(
  data = umap_scurve)

umap_scurve_scaled <- scurve_scaled_obj$scaled_nldr
lim1 <- scurve_scaled_obj$lim1
lim2 <- scurve_scaled_obj$lim2
r2 <- diff(lim2)/diff(lim1)

df_all_scurve <- dplyr::bind_cols(training_data_scurve |> dplyr::select(-ID),
                                  umap_data_with_hb_id)

# ### Define type column
# df <- df_all_scurve |>
#   dplyr::select(tidyselect::starts_with("x")) |>
#   dplyr::mutate(type = "data") ## original dataset
#
# df_b <- df_bin_scurve |>
#   dplyr::filter(hb_id %in% df_bin_centroids_scurve$hexID) |>
#   dplyr::mutate(type = "model") ## Data with summarized mean
#
# ## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
# df_b <- df_b[match(df_bin_centroids_scurve$hexID, df_b$hb_id),] |>
#   select(-type)
#
# names(df_b)[2:NCOL(df_b)] <- paste0("avg_", names(df_b)[2:NCOL(df_b)])
#
# df_all_scurve <- left_join(df_all_scurve, df_b, by = "hb_id")

shared_df_scurve <- SharedData$new(df_all_scurve |> select(-ID, -hb_id))

nldr_scurve <- shared_df_scurve |>
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha=0.5, colour="#6a3d9a", size = 0.5)

nldr_scurve_plt <- ggplotly(nldr_scurve)

langevitour_output <- langevitour(training_data |> select(-type),
                                  link=shared_df_scurve)

bscols(nldr_scurve_plt, langevitour_output)
