## Draw the axes

true_model_df <- read_rds("data/s_curve/scurve_true_model.rds")


true_model_df_rm_id <- true_model_df |>
  select(-ID)

projection <- cbind(
  c(0.18390,-0.18443,-0.06860,0.00654,-0.02921,0.02313,-0.08963),
  c(-0.02902,0.03638,-0.16872,0.22418,0.01146,0.02777,0.01453))

projected_df <- as.matrix(true_model_df_rm_id) %*% projection |>
  as_tibble(.name_repair = "unique") |>
  rename(c("proj1" = "...1",
           "proj2" = "...2"))


projected_df <- as.matrix(true_model_df_rm_id) %*% projection

limits <- 2
rng <- range(projected_df)
projected_df <- projected_df/max(abs(rng))
colnames(projected_df) <- c("P1", "P2")

obs_labels <- as.character(1:nrow(true_model_df_rm_id))

projected_df$obs_labels <- obs_labels

#if (position == "center") {
  axis_scale <- limits
  #axis_pos <- 0
  axis_pos <- -0.6 * limits
#}
# else if (position == "bottomleft") {
#   axis_scale <- limits/6
#   axis_pos <- -2/3 * limits
# }

mean_cols <- colMeans(true_model_df_rm_id, na.rm = TRUE)

adj <- function(mean, x) mean * x + x * axis_scale
axes <- data.frame(x1 = mean_cols * projection[, 1],
                   y1 = mean_cols * projection[, 2],
                   x2 = adj(mean_cols, projection[, 1]),
                   y2 = adj(mean_cols, projection[, 2]))

# axes <- data.frame(x1 = mean_cols + mean_cols * axis_scale,
#                    y1 = mean_cols + mean_cols * axis_scale,
#                    x2 = mean_cols + projection[, 1] * axis_scale,
#                    y2 = mean_cols + projection[, 2]* axis_scale)

axis_labels <- colnames(true_model_df_rm_id)
rownames(axes) <- axis_labels
#theta <- seq(0, 2 * pi, length = 50)
#circle <- data.frame(c1 = adj(cos(theta)), c2 = adj(sin(theta)))


ggplot() +
  geom_point(data = projected_df_n,
             aes(
               x = proj1,
               y = proj2),
             alpha = 0.5)  +
  geom_segment(data=axes, aes(x=x1, y=y1, xend=x2, yend=y2), colour="grey70") +
  geom_text(data=axes, aes(x=x2, y=y2, label=rownames(axes)), colour="grey50") +
  geom_point(data=axes, aes(x=x1, y=y1), colour="red")
