library(readr)
library(dplyr)

set.seed(20240110)

source("quollr_code.R", local = TRUE)
source("nldr_code.R", local = TRUE)

## Import data
df_s_curve <- read_rds("data/s_curve/s_curve.rds")
training_data_s_curve <- read_rds("data/s_curve/s_curve_training.rds")

TriMAP_s_curve <- read_rds("data/s_curve/s_curve_trimap.rds")

num_bins_s_curve <- 6

shape_value_s_curve <- calculate_effective_shape_value(.data = TriMAP_s_curve,
                                                       x = TriMAP1, y = TriMAP2) ## 1.259938

## To extract bin centroids
hexbin_data_object_s_curve <- extract_hexbin_centroids(nldr_df = TriMAP_s_curve,
                                                       num_bins = num_bins_s_curve,
                                                       shape_val = shape_value_s_curve, x = TriMAP1, y = TriMAP2)

df_bin_centroids_s_curve <- hexbin_data_object_s_curve$hexdf_data

TriMAP_s_curve_with_hb_id <- TriMAP_s_curve |>
  dplyr::mutate(hb_id = hexbin_data_object_s_curve$hb_data@cID)

## To generate a data set with high-D and 2D training data
df_all_s_curve <- dplyr::bind_cols(training_data_s_curve |> dplyr::select(-ID), TriMAP_s_curve_with_hb_id)

## Averaged on high-D
df_bin_s_curve <- avg_highD_data(.data = df_all_s_curve, column_start_text = "x")

## Triangulate bin centroids
tr1_object_s_curve <- triangulate_bin_centroids(df_bin_centroids_s_curve, x, y)
tr_from_to_df_s_curve <- generate_edge_info(triangular_object = tr1_object_s_curve)

## Compute 2D distances
distance_s_curve <- cal_2D_dist(.data = tr_from_to_df_s_curve)

## To find the benchmark value
benchmark <- find_benchmark_value(.data = distance_s_curve, distance_col = distance)
#benchmark <- 4.8


distance <- distance_s_curve

###############

distance <- distance %>%
  rename(distance_2D = distance)

#num_pca <- 4

#high_d <- training_data_s_curve
high_d <- df_bin_s_curve |> select(-hb_id)
#rownames(high_d) <- df_bin_s_curve$hb_id
#high_d <- tibble::rowid_to_column(high_d)

data_matrix <- data.matrix(high_d, rownames.force = NA)

# calculate Euclidean distance between each row in matrix
d <- dist(data_matrix)


df <- reshape2::melt(as.matrix(d), varnames = c("Start_point", "End_point"))


df <- df[df$Start_point > df$End_point,]

#rownames(df) <- 1:nrow(df)
df <- tibble::rowid_to_column(df)
df <- tibble::column_to_rownames(df, var = "rowid")


df <- df %>%
  rename(distance_high_D = value)

high_low_dist <- left_join(distance, df, by = c("from" = "End_point", "to" = "Start_point"))
high_low_dist <- high_low_dist %>% distinct()
high_low_dist <- high_low_dist[complete.cases(high_low_dist),]

lm_fit <- lm(distance_high_D ~ distance_2D, data = high_low_dist) #Adjusted R-squared:  0.7887

high_low_dist <- high_low_dist %>%
  mutate(Edge_type = ifelse(distance_2D >= benchmark, "Long edge", "Short edge"))

diagnostic_plot2 <- ggplot(high_low_dist, aes(x = distance_2D, y = distance_high_D, color = Edge_type)) +
  geom_point(alpha = 0.5) +
  geom_line(data = broom::augment(lm_fit),
            aes(x = distance_2D, y = .fitted),
            color="blue",
            linewidth=1) +
  geom_rug(col="#e31a1c",alpha=0.1, linewidth=1.5) +
  #ggtitle("Visualization of 2D and high-D distances \n by edge type.") +
  scale_color_manual(values=c("#33a02c","#000000")) +
  theme(legend.position = "bottom", legend.title=element_blank(), plot.title=element_text(size=8), axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        legend.text = element_text(size=7)) +
  xlab(expression(d^{(2)})) +
  ylab(expression(d^{(p)}))


trimesh_gr_s_curve_trimap <- colour_long_edges(.data = distance_s_curve, benchmark_value = benchmark,
                                              triangular_object = tr1_object_s_curve, distance_col = distance)

tour1_s_curve_trimap <- show_langevitour(df_all_s_curve, df_bin_s_curve,
                                        df_bin_centroids_s_curve, benchmark_value = benchmark,
                                        distance = distance_s_curve, distance_col = distance, col_start = "x")


a <- broom::augment(lm_fit) |> mutate(edge_type = if_else(distance_2D >= benchmark, "long edge", "short edge"))
a_short <- a |> filter(edge_type == "short edge")
error <- a_short$.resid/(a$.resid |> sum())

tibble::tibble(benchmark = benchmark, total_error = error |> sum())

##################

benchmark_dist_vec <- distance_s_curve$distance |> round(3) |> unique() |> sort()
#benchmark_dist_vec <- seq(min(distance_s_curve$distance |> round(3) |> unique()), max(distance_s_curve$distance |> round(3) |> unique()), 1)

vec <- stats::setNames(rep("", 2), c("benchmark_rm_lg", "total_error"))  ## Define column names

eval_data_training <- dplyr::bind_rows(vec)[0, ]
eval_data_training <- eval_data_training |>
  dplyr::mutate_if(is.character, as.numeric)

for(i in 1:length(benchmark_dist_vec)) {

  df_with_residuals <- broom::augment(lm_fit) |>
    mutate(edge_type = if_else(distance_2D >= benchmark_dist_vec[i], "long edge", "short edge"))
  df_with_residuals_short <- df_with_residuals |>
    filter(edge_type == "short edge")
  #tot <- (df_with_residuals_short$.resid |> sum())/(df_with_residuals$.resid |> sum())
  tot <- (df_with_residuals_short$.resid)/((df_with_residuals$.resid |> sum()) * NROW(df_with_residuals_short))

  df_res_benchamrk <- tibble::tibble(benchmark_rm_lg = benchmark_dist_vec[i], total_error = tot |> sum())

  eval_data_training <- dplyr::bind_rows(eval_data_training, df_res_benchamrk)



}

ggplot(eval_data_training, aes(x = benchmark_rm_lg,
                               y = total_error
)) +
  geom_point() +
  geom_line()

eval_data_training <- eval_data_training |>
  dplyr::mutate(method = "TriMAP")

write_rds(eval_data_training, "data/s_curve/s_curve_summary_lg_trimap.rds")
