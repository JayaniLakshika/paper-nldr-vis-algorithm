### Points only
X <- df |> dplyr::select(-type)

projection <- cbind(
  c(-0.04203,0.04487,-0.09140,0.05530,0.24361,-0.07504,-0.04631),
  c(0.00452,0.11483,-0.11715,0.18885,-0.11697,0.01529,-0.07619))
projected <- as.matrix(X) %*% projection
a <- projected |>
  as.data.frame() |>
  mutate(type = df$type)

ggplot(data = a, aes(x = V1, y = V2, color = type)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("#6a3d9a", "#33a02c", "#969696"))


### Points + lines
a <- a |>
  mutate(ID = row_number())

proj_true_model <- a |>
  filter(type == "true model")

proj_model <- a |>
  filter(type == "model")

connections_df <- left_join(connections, proj_true_model, by = c("from" = "ID"))
names(connections_df)[3:NCOL(connections_df)] <- paste0(names(proj_true_model)[-NCOL(proj_true_model)], "_from")

connections_df <- left_join(connections_df, proj_true_model, by = c("to" = "ID"))
names(connections_df)[(NCOL(connections) + NCOL(proj_true_model)):NCOL(connections_df)] <- paste0(names(proj_true_model)[-NCOL(proj_true_model)], "_to")


df_b <- df_b |>
  select(-type) |>
  mutate(ID = row_number() + max(connections$from) + 1)

model_df <- left_join(distance_df_small_edges, proj_model, by = c("from" = "ID"))
names(model_df)[3:NCOL(model_df)] <- paste0(names(proj_model)[-NCOL(proj_model)], "_from")

model_df <- left_join(model_df, proj_model, by = c("to" = "ID"))
names(model_df)[(NCOL(distance_df_small_edges) + NCOL(proj_model)):NCOL(model_df)] <- paste0(names(proj_model)[-NCOL(proj_model)], "_to")

line_df <- bind_rows(connections_df, model_df)

ggplot(data = a, aes(x = V1, y = V2, color = type)) +
  geom_point(alpha = 0.5) +
  geom_segment(data = line_df, aes(x = V1_from, y = V2_from, xend = V1_to, yend = V2_to, color = type_to)) +
  scale_color_manual(values = c("#6a3d9a", "#33a02c", "#969696"))

