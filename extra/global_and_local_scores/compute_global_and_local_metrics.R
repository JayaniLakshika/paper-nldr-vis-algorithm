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


# --- KNN Accuracy with LOOCV ---
knn_accuracy <- function(X_low, labels, k_values = 1:20) {
  # X_low: n x d matrix/data.frame of low-dimensional embedding
  # labels: factor or vector of class labels
  # k_values: range of k to tune

  # Ensure labels are factors
  labels <- as.factor(labels)

  n <- nrow(X_low)
  acc_results <- numeric(length(k_values))

  for (i in seq_along(k_values)) {
    k <- k_values[i]

    # Leave-one-out cross-validation
    preds <- sapply(1:n, function(j) {
      train_idx <- setdiff(1:n, j)
      test_idx <- j

      # Fit KNN on training set
      pred <- class::knn(train = X_low[train_idx, , drop = FALSE],
                         test  = X_low[test_idx, , drop = FALSE],
                         cl    = labels[train_idx],
                         k     = k)
      as.character(pred)
    })

    acc_results[i] <- mean(preds == as.character(labels))
  }

  # Pick best k
  best_idx <- which.max(acc_results)

  return(list(
    best_k = k_values[best_idx],
    best_accuracy = acc_results[best_idx],
    accuracy_by_k = data.frame(k = k_values, accuracy = acc_results)
  ))
}

set.seed(123)

# Fake embedding: 3 clusters in 2D
X_low <- rbind(
  matrix(rnorm(50*2, mean=0), ncol=2),
  matrix(rnorm(50*2, mean=3), ncol=2),
  matrix(rnorm(50*2, mean=-3), ncol=2)
)
labels <- rep(1:3, each=50)

# Run KNN Accuracy
res <- knn_accuracy(X_low, labels, k_values = 1:10)

print(res$best_k)
print(res$best_accuracy)
print(res$accuracy_by_k)

library(kernlab)
library(caret)

# --- SVM Accuracy with RBF kernel and 5-fold CV ---
svm_accuracy <- function(X_low, labels, folds = 5, sigma = 1, C = 1) {
  # X_low: n x d matrix or data.frame of low-dimensional embedding
  # labels: factor or vector of class labels
  # folds: number of CV folds
  # sigma: RBF kernel parameter (gamma = 1/(2*sigma^2))
  # C: regularization parameter

  labels <- as.factor(labels)
  n <- nrow(X_low)

  # Define CV folds
  set.seed(123)
  fold_ids <- sample(rep(1:folds, length.out = n))

  acc <- numeric(folds)

  for (f in 1:folds) {
    train_idx <- which(fold_ids != f)
    test_idx  <- which(fold_ids == f)

    # Train SVM with RBF kernel
    svm_model <- ksvm(as.matrix(X_low[train_idx, ]),
                      labels[train_idx],
                      kernel = "rbfdot",
                      kpar = list(sigma = sigma),
                      C = C,
                      scaled = FALSE)

    preds <- predict(svm_model, as.matrix(X_low[test_idx, ]))

    acc[f] <- mean(preds == labels[test_idx])
  }

  return(list(
    mean_accuracy = mean(acc),
    sd_accuracy   = sd(acc),
    fold_accuracies = acc
  ))
}

set.seed(42)

# Fake data: 3 clusters in 2D
X_low <- rbind(
  matrix(rnorm(50*2, mean=0), ncol=2),
  matrix(rnorm(50*2, mean=3), ncol=2),
  matrix(rnorm(50*2, mean=-3), ncol=2)
)
labels <- rep(1:3, each=50)

# Run SVM accuracy
res <- svm_accuracy(X_low, labels, folds = 5, sigma = 0.5, C = 1)

print(res$mean_accuracy)
print(res$sd_accuracy)
print(res$fold_accuracies)


