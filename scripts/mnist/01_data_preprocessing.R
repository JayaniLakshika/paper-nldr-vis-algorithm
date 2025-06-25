library(readr)
library(reticulate)
library(dplyr)
library(ggplot2)

mnist_data <- read_rds("~/Desktop/PhD Monash research files/paper-nldr-vis-algorithm/data/mnist/mnist_digit_1.rds")

data_x <- mnist_data[,-785]/255
data_x <- array_reshape(data_x, c(nrow(data_x), 784), order = "F")
data_x_df <- as.data.frame(data_x)
names(data_x_df) <- paste0(rep("x", ncol(data_x_df)), 1:ncol(data_x_df))

# PCA

# data_x_df <- data_x_df |>
#   select(where(~ any(. != 0)))

dim(data_x_df)


#Spectral decomposition which examines the covariances / correlations between variables
#princomp(data, cor = FALSE)
calculate_pca <- function(feature_dataset, num_pcs){
  pcaY_cal <- prcomp(feature_dataset, center = TRUE, scale. = FALSE)
  PCAresults <- data.frame(pcaY_cal$x[, 1:num_pcs])
  summary_pca <- summary(pcaY_cal)
  var_explained_df <- data.frame(PC= paste0("PC",1:50),
                                 var_explained=(pcaY_cal$sdev[1:50])^2/sum((pcaY_cal$sdev[1:50])^2))
  return(list(prcomp_out = pcaY_cal,pca_components = PCAresults, summary = summary_pca, var_explained_pca  = var_explained_df))
}

pca_ref_calc <- calculate_pca(data_x_df,10)
pca_ref_calc$summary
var_explained_df <- pca_ref_calc$var_explained_pca
data_pca <- pca_ref_calc$pca_components

var_explained_df |>
  ggplot(aes(x = PC,y = var_explained, group = 1))+
  geom_point(size=1)+
  geom_line()+
  labs(title="Scree plot: PCA on scaled data") +
  scale_x_discrete(limits = paste0(rep("PC", 50), 1:50)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

write_rds(data_pca, "data/mnist/mnist_10_pcs_of_digit_1.rds")


