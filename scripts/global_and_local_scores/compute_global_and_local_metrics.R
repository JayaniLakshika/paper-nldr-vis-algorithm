## high-dimensional data (X_high), its embedding (X_low), and optionally labels (y).

# --- Utility: Triplet Accuracy for a given set of triplets ---
triplet_accuracy <- function(X_high, X_low, triplets) {
  n_triplets <- nrow(triplets)

  # Compute distances in high- and low-dimensional space
  d_high <- function(i, j) sum((X_high[i, ] - X_high[j, ])^2)
  d_low  <- function(i, j) sum((X_low[i, ] - X_low[j, ])^2)

  correct <- 0
  for (t in 1:n_triplets) {
    i <- triplets[t, 1]
    j <- triplets[t, 2]
    k <- triplets[t, 3]

    # Order in high-D
    d_ij_high <- d_high(i, j)
    d_ik_high <- d_high(i, k)
    high_closer <- ifelse(d_ij_high < d_ik_high, "j", "k")

    # Order in low-D
    d_ij_low <- d_low(i, j)
    d_ik_low <- d_low(i, k)
    low_closer <- ifelse(d_ij_low < d_ik_low, "j", "k")

    if (high_closer == low_closer) {
      correct <- correct + 1
    }
  }
  return(correct / n_triplets)
}


# --- Random Triplet Accuracy (RTA) ---
random_triplet_accuracy <- function(X_high, X_low, n_triplets = 10000, n_repeats = 5) {
  n <- nrow(X_high)
  results <- numeric(n_repeats)

  for (r in 1:n_repeats) {
    # Sample triplets (i, j, k)
    triplets <- t(replicate(n_triplets, sample(1:n, 3, replace = FALSE)))
    acc <- triplet_accuracy(X_high, X_low, triplets)
    results[r] <- acc
  }

  return(list(mean = mean(results), sd = sd(results), values = results))
}


# --- Centroid Triplet Accuracy (CTA) ---
centroid_triplet_accuracy <- function(X_high, X_low, labels) {
  # Compute centroids
  centroids_high <- aggregate(X_high, by = list(labels), FUN = mean)
  centroids_low  <- aggregate(X_low,  by = list(labels), FUN = mean)

  # Drop label column
  centroids_high <- as.matrix(centroids_high[, -1])
  centroids_low  <- as.matrix(centroids_low[, -1])

  n <- nrow(centroids_high)
  triplets <- t(combn(n, 3)) # all possible centroid triplets

  acc <- triplet_accuracy(centroids_high, centroids_low, triplets)
  return(acc)
}

set.seed(123)

# Fake data: high-D (5D), low-D (2D), 3 classes
X_high <- matrix(rnorm(300 * 5), ncol = 5)
X_low <- matrix(rnorm(300 * 2), ncol = 2)
labels <- sample(1:3, 300, replace = TRUE)

# Random Triplet Accuracy
rta <- random_triplet_accuracy(X_high, X_low, n_triplets = 5000, n_repeats = 5)
print(rta)
# $mean and $sd will give summary

# Centroid Triplet Accuracy
cta <- centroid_triplet_accuracy(X_high, X_low, labels)
print(cta)

