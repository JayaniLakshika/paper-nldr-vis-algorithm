# Set the number of samples
n_samples <- 5000  # Replace with your desired number of samples
noise <- 0.0      # Replace with your desired noise level

# Generate the t values
t <- 3 * pi * (runif(n_samples) - 0.5)

# Create an empty matrix to store the results
X <- matrix(nrow = n_samples, ncol = 3)

# Fill in the X matrix with the corresponding values
X[, 1] <- sin(t)
X[, 2] <- 2.0 * runif(n_samples)
X[, 3] <- sign(t) * (cos(t) - 1)

# Add noise
X <- X + noise * matrix(rnorm(3 * n_samples), nrow = n_samples, ncol = 3)

# If needed, t can be squeezed (not necessary in R since t is already a vector)
t <- t

df <- as_tibble(X)
names(df) <- paste0("x", 1:3)

df$x4 <- runif(n_samples, -0.02, 0.02)
df$x5 <- runif(n_samples, -0.02, 0.02)
df$x6 <- runif(n_samples, -0.01, 0.01)
df$x7 <- runif(n_samples, -0.01, 0.01)
langevitour(df)

training_data <- df |>
  dplyr::mutate(type = "data")
