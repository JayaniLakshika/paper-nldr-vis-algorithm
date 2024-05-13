### Define type column
df <- df_all_pbmc |>
  dplyr::select(tidyselect::starts_with("PC")) |>
  dplyr::mutate(type = "data") ## original dataset

df_b <- df_bin_pbmc |>
  dplyr::filter(hb_id %in% df_bin_centroids_pbmc$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_pbmc$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

df_exe <- dplyr::bind_rows(df_b, df)

## Set the maximum difference as the criteria
distance_df_small_edges <- distance_pbmc |>
  dplyr::filter(distance < benchmark_pbmc)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to,
                         group = df_exe$type, pointSize = append(rep(0, NROW(df_b)), rep(0.5, NROW(df))),
                         levelColors = c("#6a3d9a", "#33a02c"))
