mnist_data <- read_rds("data/mnist/mnist_digit_1.rds")

img_right_top <- c(6561, 454, 6397)

pixels_gathered_within <-  mnist_data |>
  mutate(instance = row_number()) |>
  gather(pixel, value, -Label, -instance) |>
  extract(pixel, "pixel", "(\\d+)", convert = TRUE) |>
  mutate(pixel = pixel - 2, x = pixel %% 28, y = 28 - pixel %/% 28) |>
  filter(instance %in% img_error_within)

right_top_img <- pixels_gathered_within |>
  ggplot(aes(x, y, fill = value)) +
  geom_tile() +
  facet_wrap(~ instance, ncol = 4) +
  coord_fixed() +
  scale_fill_continuous_sequential(palette = "Grays") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")

img_middle <- c(4671, 3444, 3694)

pixels_gathered_within <-  mnist_data |>
  mutate(instance = row_number()) |>
  gather(pixel, value, -Label, -instance) |>
  extract(pixel, "pixel", "(\\d+)", convert = TRUE) |>
  mutate(pixel = pixel - 2, x = pixel %% 28, y = 28 - pixel %/% 28) |>
  filter(instance %in% img_error_within)

middle_img <- pixels_gathered_within |>
  ggplot(aes(x, y, fill = value)) +
  geom_tile() +
  facet_wrap(~ instance, ncol = 4) +
  coord_fixed() +
  scale_fill_continuous_sequential(palette = "Grays") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")

img_right_bottom <- c(296, 6386, 2)

pixels_gathered_within <-  mnist_data |>
  mutate(instance = row_number()) |>
  gather(pixel, value, -Label, -instance) |>
  extract(pixel, "pixel", "(\\d+)", convert = TRUE) |>
  mutate(pixel = pixel - 2, x = pixel %% 28, y = 28 - pixel %/% 28) |>
  filter(instance %in% img_error_within)

right_bottom_img <- pixels_gathered_within |>
  ggplot(aes(x, y, fill = value)) +
  geom_tile() +
  facet_wrap(~ instance, ncol = 4) +
  coord_fixed() +
  scale_fill_continuous_sequential(palette = "Grays") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")
