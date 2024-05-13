### Define type column
df <- df_all_mnist |>
  dplyr::select(tidyselect::starts_with("PC")) |>
  dplyr::mutate(type = "data") ## original dataset

df_b <- df_bin_mnist |>
  dplyr::filter(hb_id %in% df_bin_centroids_mnist$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_mnist$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

df_exe <- dplyr::bind_rows(df_b, df)

## Set the maximum difference as the criteria
distance_df_small_edges <- distance_mnist |>
  dplyr::filter(distance < benchmark_mnist)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to,
                         group = df_exe$type, pointSize = append(rep(0, NROW(df_b)), rep(0.5, NROW(df))),
                         levelColors = c("#6a3d9a", "#33a02c"))

### With error
error_df_n <- error_df |>
  select(starts_with("PC"), type)

df_b <- df_bin_mnist |>
  dplyr::filter(hb_id %in% df_bin_centroids_mnist$hexID) |>
  dplyr::mutate(type = "model") ## Data with summarized mean

## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
df_b <- df_b[match(df_bin_centroids_mnist$hexID, df_b$hb_id),] |>
  dplyr::select(-hb_id)

df_exe <- dplyr::bind_rows(df_b, error_df_n)

## Set the maximum difference as the criteria
distance_df_small_edges <- distance_mnist |>
  dplyr::filter(distance < benchmark_mnist)
## Since erase brushing is considerd.

langevitour::langevitour(df_exe[1:(length(df_exe)-1)],
                         lineFrom = distance_df_small_edges$from,
                         lineTo = distance_df_small_edges$to,
                         group = factor(df_exe$type, levels = c("error 0-2",
                                                                "error 2-4", "error 4-6",
                                                                "error 6-8",
                                                                "error 8-10", "error 10-12",
                                                                "error 12-14", "error greter than 14",
                                                                "model")),
                         pointSize = append(rep(0, NROW(df_b)), rep(1, NROW(df))),
                         levelColors = append(rainbow_hcl(8), "#33a02c"))

