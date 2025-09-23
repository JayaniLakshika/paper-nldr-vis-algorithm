## This script is to compute different evaluation metrics suggested by PaCMAP and TriMAP authors

# Load the reticulate library
library(reticulate)
use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

# Run the Python script
py_run_file("~/Desktop/PhD Monash research files/Research papers/paper-nldr-vis-algorithm/scripts/evaluation.py")

######Global score with TriMAP #############
# Create an instance of the Python class
my_instance <- py$MyClass()

# Assuming you have R matrices for X and Y
X_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
Y_matrix <- matrix(rnorm(50), nrow = 10, ncol = 5)

# Pass R matrices to the Python method; they will be converted to NumPy arrays
score_result <- my_instance$global_score(X_matrix, Y_matrix)
print(score_result)

###### PaCMAP evaluations ##############

# Create sample data in R
X_data <- matrix(rnorm(350), nrow = 35, ncol = 10)
X_new_data <- matrix(rnorm(70), nrow = 35, ncol = 2)
y_data <- sample(0:1, 35, replace = TRUE)

py$random_triplet_eval(X = X_data, X_new = X_new_data, y = y_data) #working

# Call the Python function and pass the R data
# The arguments correspond to the function's signature: X, X_new, y, name, baseline, labelled
results <- py$evaluate_output(X = X_data, X_new = X_new_data, y = y_data, name = "my_run", baseline = FALSE, labelled = TRUE)

# The result is an R list
print(results)
