## Data with pixel values
mnist_data <- read_rds("data/mnist/mnist_digit_1.rds")

pixels_gathered <-  mnist_data |>
  mutate(instance = row_number()) |>
  gather(pixel, value, -Label, -instance) |>
  tidyr::extract(pixel, "pixel", "(\\d+)", convert = TRUE) |>
  mutate(pixel = pixel - 2, x = pixel %% 28, y = 28 - pixel %/% 28) |>
  filter(instance %in% c(3710, 2391, 2385, 5030, 6475,
                         7679, 3568, 1312, 3311, 1097,
                         3552, 7853, 6489, 7689, 6690,
                         1380, 6057, 2347, 5946, 3355,
                         4175, 3997, 5378, 387, 1854,
                         614, 3079, 1762, 5239, 3723,
                         5748, 728, 7419, 7794, 6233,
                         33, 2485, 5998, 318, 1761,
                         5690, 165, 517, 6935, 2682,
                         4962, 2264, 5563, 6369, 559,
                         5188, 4849, 2666, 3448, 3055,
                         3120, 6869, 6345, 4470, 7147)) |>
  mutate(instance.labs =  case_when(
    instance %in% c(3710, 2391, 2385, 5030, 6475) ~ "1",
    instance %in% c(7679, 3568, 1312, 3311, 1097) ~ "2",
    instance %in% c(3552, 7853, 6489, 7689, 6690) ~ "3",
    instance %in% c(1380, 6057, 2347, 5946, 3355) ~ "4",
    instance %in% c(4175, 3997, 5378, 387, 1854) ~ "5",
    instance %in% c(614, 3079, 1762, 5239, 372) ~ "6",
    instance %in% c(5748, 728, 7419, 7794, 6233) ~ "7",
    instance %in% c(5690, 165, 517, 6935, 2682) ~ "8",
    instance %in% c(4962, 2264, 5563, 6369, 559) ~ "9",
    instance %in% c(5188, 4849, 2666, 3448, 3055) ~ "10",
    instance %in% c(3120, 6869, 6345, 4470, 7147) ~ "11",
    .default = "12"
  ))

## To generate labels
instance.labs <- as.character(rep(1:12, each = 5))

img_sample <- pixels_gathered |>
  ggplot(aes(x, y, fill = value)) +
  geom_tile() +
  facet_wrap(~ factor(instance, levels = c(3710, 2391, 2385, 5030, 6475,
                                           7679, 3568, 1312, 3311, 1097,
                                           3552, 7853, 6489, 7689, 6690,
                                           1380, 6057, 2347, 5946, 3355,
                                           4175, 3997, 5378, 387, 1854,
                                           614, 3079, 1762, 5239, 3723,
                                           5748, 728, 7419, 7794, 6233,
                                           33, 2485, 5998, 318, 1761,
                                           5690, 165, 517, 6935, 2682,
                                           4962, 2264, 5563, 6369, 559,
                                           5188, 4849, 2666, 3448, 3055,
                                           3120, 6869, 6345, 4470, 7147)), nrow = 5) +
  coord_fixed() +
  scale_fill_continuous_sequential(palette = "Grays") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")
