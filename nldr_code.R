calculate_effective_perplexity <- function(data){

  opt_perplexity <- data %>%
    NROW() %>%
    sqrt()

  round(opt_perplexity,0)

}

# num_dim: Number of dimensions to fit (default is 2): 2, or 3

Fit_tSNE <- function(data, opt_perplexity = NA, num_dim = 2, with_seed = NULL, ...){
  # If a seed is specified, then use it, otherwise ignore
  if(!is.null(with_seed)){set.seed(with_seed)}

  if(is.na(opt_perplexity)){
    opt_perplexity <- calculate_effective_perplexity(data)
  }

  tSNE_fit <- data %>%
    select(where(is.numeric)) %>%
    Rtsne::Rtsne(perplexity = opt_perplexity, pca = FALSE, pca_center = FALSE, normalize = FALSE, dims = num_dim)

  tSNE_df <- tSNE_fit$Y %>%
    tibble::as_tibble(.name_repair = "unique")  %>%
    mutate(ID = row_number())

  names(tSNE_df)[1:(NCOL(tSNE_df) - 1)] <- paste0(rep("tSNE",(NCOL(tSNE_df) - 1)), 1:(NCOL(tSNE_df) - 1))
  return(tSNE_df)
}

plot_tSNE_2D <- function(tSNE_df){
  tSNE_df_plot <- tSNE_df %>%
    ggplot(aes(x = tSNE1,
               y = tSNE2))+
    geom_point(alpha=0.5, colour="#e41a1c", size = 0.5) +
    coord_equal() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text = element_text(size = 5),
          axis.title = element_text(size = 7))
  return(tSNE_df_plot)
}


# num_dim: Number of dimensions to fit (default is 2): 2, or 3

Fit_UMAP <- function(data, n_neighbors = NA, num_dim = 2, with_seed = NULL, ...){
  # If a seed is specified, then use it, otherwise ignore
  if(!is.null(with_seed)){set.seed(with_seed)}

  if(is.na(n_neighbors)){
    n_neighbors <- 15
  }

  UMAP_fit <- data %>%
    select(where(is.numeric)) %>%
    umap(n_neighbors = n_neighbors, n_components =  num_dim)

  UMAP_df <- UMAP_fit$layout %>%
    as.data.frame()  %>%
    mutate(ID=row_number())

  names(UMAP_df)[1:(ncol(UMAP_df)-1)] <- paste0(rep("UMAP",(ncol(UMAP_df)-1)), 1:(ncol(UMAP_df)-1))
  return(list(UMAP_fit = UMAP_fit, UMAP_df = UMAP_df))
}

plot_UMAP_2D <- function(UMAP_df){
  UMAP_df_plot <- UMAP_df %>%
    ggplot(aes(x = UMAP1,
               y = UMAP2))+
    geom_point(alpha=0.5, colour="#377eb8", size = 0.5) +
    coord_equal() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text = element_text(size = 5),
          axis.title = element_text(size = 7))
  return(UMAP_df_plot)
}

Fit_PHATE <- function(data, knn = NA, with_seed = NULL, ...){
  # If a seed is specified, then use it, otherwise ignore
  if(!is.null(with_seed)){set.seed(with_seed)}

  if(is.na(knn)){
    knn <- 5
  }

  tree_phate_fit <- phate(data, knn)

  PHATE_df <- as.data.frame(tree_phate_fit$embedding) %>%
    mutate(ID=row_number())

  names(PHATE_df)[1:(ncol(PHATE_df)-1)] <- paste0(rep("PHATE",(ncol(PHATE_df)-1)), 1:(ncol(PHATE_df)-1))
  PHATE_df
}

plot_PHATE_2D <- function(PHATE_df){
  PHATE_df_plot <- PHATE_df %>%
    ggplot(aes(x = PHATE1,
               y = PHATE2))+
    geom_point(alpha=0.5, colour="#4daf4a", size = 0.5) +
    coord_equal() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text = element_text(size = 5),
          axis.title = element_text(size = 7))
  return(PHATE_df_plot)
}

Fit_TriMAP_data <- function(df_2_without_class, tem_dir){
  write.csv(df_2_without_class, file.path(tem_dir, "df_2_without_class.csv"), row.names = FALSE,
            quote = TRUE)
}

plot_TriMAP_2D <- function(TRIMAP_df){
  TriMAP_df_plot <- TRIMAP_df %>%
    ggplot(aes(x = TriMAP1,
               y = TriMAP2))+
    geom_point(alpha=0.5, colour="#984ea3", size = 0.5) +
    coord_equal() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text = element_text(size = 5),
          axis.title = element_text(size = 7))
  return(TriMAP_df_plot)
}

Fit_PacMAP_data <- function(df_2_without_class, tem_dir){
  write.csv(df_2_without_class, file.path(tem_dir, "df_2_without_class.csv"), row.names = FALSE,
            quote = TRUE)
}

plot_PaCMAP_2D <- function(PaCMAP_df){
  PaCMAP_df_plot <- PaCMAP_df %>%
    ggplot(aes(x = PaCMAP1,
               y = PaCMAP2))+
    geom_point(alpha=0.5, colour="#ff7f00", size = 0.5) +
    coord_equal() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text = element_text(size = 5),
          axis.title = element_text(size = 7))
  return(PaCMAP_df_plot)
}
