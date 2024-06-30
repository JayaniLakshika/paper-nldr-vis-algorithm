
<!-- five Gaussian clusters with different NLDR techniques-->
  ```{r}
#| warning: false
#| echo: false

## Import data
training_data_gau <- read_rds("data/five_gau_clusters/data_five_gau_training.rds")

tsne_data_gau <- read_rds("data/five_gau_clusters/tsne_data_five_gau_61.rds")
umap_data_gau <- read_rds("data/five_gau_clusters/umap_data_five_gau.rds")
phate_data_gau <- read_rds("data/five_gau_clusters/phate_data_five_gau.rds")
trimap_data_gau <- read_rds("data/five_gau_clusters/trimap_data_five_gau.rds")
pacmap_data_gau <- read_rds("data/five_gau_clusters/pacmap_data_five_gau.rds")

## Visualize embedding

nldr_gau1 <- tsne_data_gau |>
  ggplot(aes(x = tSNE1,
             y = tSNE2)) +
  geom_point(alpha=0.1, size=1) +
  interior_annotation("a")

nldr_gau2 <- umap_data_gau |>
  ggplot(aes(x = UMAP1,
             y = UMAP2)) +
  geom_point(alpha=0.1, size=1) +
  interior_annotation("b")

nldr_gau3 <- phate_data_gau |>
  ggplot(aes(x = PHATE1,
             y = PHATE2)) +
  geom_point(alpha=0.1, size=1) +
  interior_annotation("c")

nldr_gau4 <- pacmap_data_gau |>
  ggplot(aes(x = PaCMAP1,
             y = PaCMAP2)) +
  geom_point(alpha=0.1, size=1) +
  interior_annotation("d")

nldr_gau5 <- trimap_data_gau |>
  ggplot(aes(x = TriMAP1,
             y = TriMAP2)) +
  geom_point(alpha=0.1, size=1) +
  interior_annotation("e")
```

```{r}
#| label: fig-NLDR-variety-gau
#| echo: false
#| fig-cap: "Five different NLDR representations of the same data. Different techniques and different parameter choices are used. Is there a best representation of the original data or are they all providing  equivalent information?"
#| fig-width: 8
#| fig-height: 2
#| out-width: 100%
#| fig-pos: H
# (a) tSNE (perplexity = 61), (b) UMAP (n_neighbors = 15), (c) PHATE (knn = 5), (d) TriMAP (n_inliers = 5, n_outliers = 4, n_random = 3), and (e) PaCMAP (n_neighbors = 10, init = random, MN_ratio = 0.9, FP_ratio = 2)
nldr_gau1 + nldr_gau2 +
  nldr_gau3 + nldr_gau4 + nldr_gau5 +
  plot_layout(ncol = 5)
```

```{r}
sc_xlims_gau <- c(-0.25, 1.25)
sc_ylims_gau <- c(-0.2, 1.2)
```

<!--trimesh with tsne-->
  ```{r}
gau1_scaled_obj <- gen_scaled_data(
  data = tsne_data_gau)
tsne_gau_scaled <- gau1_scaled_obj$scaled_nldr

## Compute hexbin parameters
num_bins_x_gau1 <- 12
lim1 <- gau1_scaled_obj$lim1
lim2 <- gau1_scaled_obj$lim2
r2_gau1 <- diff(lim2)/diff(lim1)

gau1_model <- fit_highd_model(
  training_data = training_data_gau,
  emb_df = tsne_gau_scaled,
  bin1 = num_bins_x_gau1,
  r2 = r2_gau1,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_gau1 <- gau1_model$df_bin_centroids
df_bin_gau1 <- gau1_model$df_bin

## Triangulate bin centroids
tr1_object_gau1 <- tri_bin_centroids(
  df_bin_centroids_gau1, x = "c_x", y = "c_y")
tr_from_to_df_gau1 <- gen_edges(
  tri_object = tr1_object_gau1)

## Compute 2D distances
distance_gau1 <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_gau1,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_gau1 <- find_lg_benchmark(
  distance_edges = distance_gau1,
  distance_col = "distance")

trimesh_removed_gau1 <- vis_rmlg_mesh(
  distance_edges = distance_gau1,
  benchmark_value = benchmark_gau1,
  tr_coord_df = tr_from_to_df_gau1,
  distance_col = "distance") +
  xlim(sc_xlims_gau) + ylim(sc_ylims_gau) +
  interior_annotation("a", sc_ltr_pos)
```

<!--trimesh with umap-->
  ```{r}
gau2_scaled_obj <- gen_scaled_data(
  data = umap_data_gau)
umap_gau_scaled <- gau2_scaled_obj$scaled_nldr

## Compute hexbin parameters
num_bins_x_gau2 <- 44
lim1 <- gau2_scaled_obj$lim1
lim2 <- gau2_scaled_obj$lim2
r2_gau2 <- diff(lim2)/diff(lim1)

gau2_model <- fit_highd_model(
  training_data = training_data_gau,
  emb_df = umap_gau_scaled,
  bin1 = num_bins_x_gau2,
  r2 = r2_gau2,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_gau2 <- gau2_model$df_bin_centroids
df_bin_gau2 <- gau2_model$df_bin

## Triangulate bin centroids
tr1_object_gau2 <- tri_bin_centroids(
  df_bin_centroids_gau2, x = "c_x", y = "c_y")
tr_from_to_df_gau2 <- gen_edges(
  tri_object = tr1_object_gau2)

## Compute 2D distances
distance_gau2 <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_gau2,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_gau2 <- find_lg_benchmark(
  distance_edges = distance_gau2,
  distance_col = "distance")

trimesh_removed_gau2 <- vis_rmlg_mesh(
  distance_edges = distance_gau2,
  benchmark_value = benchmark_gau2,
  tr_coord_df = tr_from_to_df_gau2,
  distance_col = "distance") +
  xlim(sc_xlims_gau) + ylim(sc_ylims_gau) +
  interior_annotation("b", sc_ltr_pos)
```

<!--trimesh with phate-->
  ```{r}
# gau3_scaled_obj <- gen_scaled_data(
#   data = phate_data_gau)
# phate_gau_scaled <- gau3_scaled_obj$scaled_nldr
#
# ## Compute hexbin parameters
# num_bins_x_gau3 <- 250
# lim1 <- gau3_scaled_obj$lim1
# lim2 <- gau3_scaled_obj$lim2
# r2_gau3 <- diff(lim2)/diff(lim1)
#
# gau3_model <- fit_highd_model(
#   training_data = training_data_gau,
#   emb_df = phate_gau_scaled,
#   bin1 = num_bins_x_gau3,
#   r2 = r2_gau3,
#   is_bin_centroid = TRUE,
#   is_rm_lwd_hex = FALSE,
#   col_start_highd = "x"
# )
#
# df_bin_centroids_gau3 <- gau3_model$df_bin_centroids
# df_bin_gau3 <- gau3_model$df_bin

df_bin_centroids_gau3 <- read_rds("data/five_gau_clusters/df_bin_centroids_gau_phate.rds")
df_bin_gau3 <- read_rds("data/five_gau_clusters/df_bin_gau_phate.rds")

## Triangulate bin centroids
tr1_object_gau3 <- tri_bin_centroids(
  df_bin_centroids_gau3, x = "c_x", y = "c_y")
tr_from_to_df_gau3 <- gen_edges(
  tri_object = tr1_object_gau3)

## Compute 2D distances
distance_gau3 <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_gau3,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_gau3 <- find_lg_benchmark(
  distance_edges = distance_gau3,
  distance_col = "distance")

trimesh_removed_gau3 <- vis_rmlg_mesh(
  distance_edges = distance_gau3,
  benchmark_value = benchmark_gau3,
  tr_coord_df = tr_from_to_df_gau3,
  distance_col = "distance") +
  xlim(sc_xlims_gau) + ylim(sc_ylims_gau) +
  interior_annotation("c", sc_ltr_pos)
```

<!--trimesh with pacmap-->
  ```{r}
gau4_scaled_obj <- gen_scaled_data(
  data = pacmap_data_gau)
pacmap_gau_scaled <- gau4_scaled_obj$scaled_nldr

## Compute hexbin parameters
num_bins_x_gau4 <- 18
lim1 <- gau4_scaled_obj$lim1
lim2 <- gau4_scaled_obj$lim2
r2_gau4 <- diff(lim2)/diff(lim1)

gau4_model <- fit_highd_model(
  training_data = training_data_gau,
  emb_df = pacmap_gau_scaled,
  bin1 = num_bins_x_gau4,
  r2 = r2_gau4,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_gau4 <- gau4_model$df_bin_centroids
df_bin_gau4 <- gau4_model$df_bin

## Triangulate bin centroids
tr1_object_gau4 <- tri_bin_centroids(
  df_bin_centroids_gau4, x = "c_x", y = "c_y")
tr_from_to_df_gau4 <- gen_edges(
  tri_object = tr1_object_gau4)

## Compute 2D distances
distance_gau4 <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_gau4,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_gau4 <- find_lg_benchmark(
  distance_edges = distance_gau4,
  distance_col = "distance")

trimesh_removed_gau4 <- vis_rmlg_mesh(
  distance_edges = distance_gau4,
  benchmark_value = benchmark_gau4,
  tr_coord_df = tr_from_to_df_gau4,
  distance_col = "distance") +
  xlim(sc_xlims_gau) + ylim(sc_ylims_gau) +
  interior_annotation("d", sc_ltr_pos)
```

<!--trimesh with trimap-->
  ```{r}
gau5_scaled_obj <- gen_scaled_data(
  data = trimap_data_gau)
trimap_gau_scaled <- gau5_scaled_obj$scaled_nldr

## Compute hexbin parameters
num_bins_x_gau5 <- 25
lim1 <- gau5_scaled_obj$lim1
lim2 <- gau5_scaled_obj$lim2
r2_gau5 <- diff(lim2)/diff(lim1)

gau5_model <- fit_highd_model(
  training_data = training_data_gau,
  emb_df = trimap_gau_scaled,
  bin1 = num_bins_x_gau5,
  r2 = r2_gau5,
  is_bin_centroid = TRUE,
  is_rm_lwd_hex = FALSE,
  col_start_highd = "x"
)

df_bin_centroids_gau5 <- gau5_model$df_bin_centroids
df_bin_gau5 <- gau5_model$df_bin

## Triangulate bin centroids
tr1_object_gau5 <- tri_bin_centroids(
  df_bin_centroids_gau5, x = "c_x", y = "c_y")
tr_from_to_df_gau5 <- gen_edges(
  tri_object = tr1_object_gau5)

## Compute 2D distances
distance_gau5 <- cal_2d_dist(
  tr_coord_df = tr_from_to_df_gau5,
  start_x = "x_from",
  start_y = "y_from",
  end_x = "x_to",
  end_y = "y_to",
  select_vars = c("from", "to", "distance"))

## To find the benchmark value
benchmark_gau5 <- find_lg_benchmark(
  distance_edges = distance_gau5,
  distance_col = "distance")
benchmark_gau5 <- 0.2

trimesh_removed_gau5 <- vis_rmlg_mesh(
  distance_edges = distance_gau5,
  benchmark_value = benchmark_gau5,
  tr_coord_df = tr_from_to_df_gau5,
  distance_col = "distance") +
  xlim(sc_xlims_gau) + ylim(sc_ylims_gau) +
  interior_annotation("e", sc_ltr_pos)
```

<!--add trimesh for all five-->

  ```{r}
#| label: fig-trimesh-gau
#| echo: false
#| fig-cap: "Model generated with five different NLDR methods in $2\\text{-}D$ with approximately $65$ non-empty bins in each."
#| fig-width: 8
#| fig-height: 2
#| out-width: 100%
#| fig-pos: H

trimesh_removed_gau1 + trimesh_removed_gau2 +
  trimesh_removed_gau3 + trimesh_removed_gau4 + trimesh_removed_gau5 +
  plot_layout(ncol = 5)
```

<!--add three screenshots from each technique with you tube links-->

  To investigate which is the reasonable representation to visualize the five spherical Gaussian cluster data or all NLDR methods provide equivalent information, we visualize all the models in \pD{} space. Models from all NLDR methods show five well-separated clusters (see @fig-gau1-sc, @fig-gau2-sc, @fig-gau3-sc, @fig-gau4-sc, and @fig-gau5-sc). This suggests that for the five Gaussian cluster dataset, all NLDR methods effectively preserve the global structure. tSNE displays clusters with varying densities, indicating their ability to capture within-cluster variation (see @fig-gau1-sc). On the other hand, both UMAP, PHATE, PaCMAP and TriMAP show clusters with flat surfaces, suggesting a failure to capture within-cluster variation (see @fig-gau2-sc, @fig-gau3-sc, @fig-gau4-sc and @fig-gau5-sc). Therefore, UMAP, PHATE, PaCMAP and TriMAP do not capture the local structure as effectively as other methods.

::: {#fig-gau1-sc layout-ncol="3" fig-pos="H"}
  ![](figures/five_gau_clusters/sc_tsne_1.png)

  ![](figures/five_gau_clusters/sc_tsne_2.png)

  ![](figures/five_gau_clusters/sc_tsne_3.png)

  Screen shots of the **langevitour** of the five Gaussian clusters dataset, shows the model with tSNE in \pD{}, a video of the tour animation is available at (<https://youtu.be/RASEE7N5MbM>).
  :::

    ::: {#fig-gau2-sc layout-ncol="3" fig-pos="H"}
      ![](figures/five_gau_clusters/sc_umap_1.png)

      ![](figures/five_gau_clusters/sc_umap_2.png)

      ![](figures/five_gau_clusters/sc_umap_3.png)

      Screen shots of the **langevitour** of the five Gaussian clusters dataset, shows the model with UMAP in \pD{}, a video of the tour animation is available at (<https://youtu.be/iG4bCPkJilw>).
      :::

        ::: {#fig-gau3-sc layout-ncol="3" fig-pos="H"}
          ![](figures/five_gau_clusters/sc_phate_1.png)

          ![](figures/five_gau_clusters/sc_phate_2.png)

          ![](figures/five_gau_clusters/sc_phate_3.png)

          Screen shots of the **langevitour** of the five Gaussian clusters dataset, shows the model with PHATE in \pD{}, a video of the tour animation is available at (<https://youtu.be/L_PVLGwfOS0>).
          :::

            ::: {#fig-gau4-sc layout-ncol="3" fig-pos="H"}
              ![](figures/five_gau_clusters/sc_pacmap_1.png)

              ![](figures/five_gau_clusters/sc_pacmap_2.png)

              ![](figures/five_gau_clusters/sc_pacmap_3.png)

              Screen shots of the **langevitour** of the five Gaussian clusters dataset, shows the model with PaCMAP in \pD{}, a video of the tour animation is available at (<https://youtu.be/z07cKXi8EJQ>).
              :::

                ::: {#fig-gau5-sc layout-ncol="3" fig-pos="H"}
                  ![](figures/five_gau_clusters/sc_trimap_1.png)

                  ![](figures/five_gau_clusters/sc_trimap_2.png)

                  ![](figures/five_gau_clusters/sc_trimap_3.png)

                  Screen shots of the **langevitour** of the five Gaussian clusters dataset, shows the model with TriMAP in \pD{}, a video of the tour animation is available at (<https://youtu.be/Chs1lYAoX2w>).
                  :::

                    When compare the NLDR representations and generated models, tSNE with perplexity $61$ appears to be a reasonable representation for visualizing the five Gaussian cluster dataset. This is supported by investigating the model generated with tSNE in the data space, which provides evidence that it preserves both local and global structures. Also, the NLDR representation with tSNE shows five well-separated clusters.
