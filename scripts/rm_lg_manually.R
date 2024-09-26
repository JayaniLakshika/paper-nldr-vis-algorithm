d <- tr_from_to_df_gau1 |>
  mutate(ID = row_number())

a <- ggplot(data = d,
            aes(
              x = x_from,
              y = y_from,
              xend = x_to,
              yend = y_to,
              group = ID)) +
  geom_segment(colour = "#000000") +
  geom_text(aes(label = from),
            size = 2)

ggplotly(a)

