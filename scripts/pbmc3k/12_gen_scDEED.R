#For scDEED method: https://github.com/JSB-UCLA/scDEED?tab=readme-ov-file
# For PBMC3k data:https://github.com/XiDsLab/Festem_paper/tree/main

library(scDEED)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)
set.seed(20240110)

## Updated version of the scDEED function to store the embedding as well (useful to compare with our method)

Distances.tSNE_add_seu_obj <- function (pbmc, pbmc.permuted, K, perplexity_score = 40, pre_embedding = "pca",
                            check_duplicates = T, rerun = T)
{
  distances <- distances::distances
  if (rerun) {
    pbmc = Seurat::RunTSNE(pbmc, seed.use = 100, perplexity = perplexity_score,
                           reduction = pre_embedding, do.fast = T, check_duplicates = check_duplicates)
  }
  ## To obtain tSNE embeddings
  tSNE_data <- pbmc@reductions$tsne@cell.embeddings
  ## write to the folder
  write_rds(tSNE_data, here::here(paste0("data/pbmc3k/pbmc_scdeed_tsne_perplexity_", perplexity_score, ".rds")))

  tSNE_distances = distances(tSNE_data)
  pbmc.permuted <- Seurat::RunTSNE(pbmc.permuted, seed.use = 100,
                                   perplexity = perplexity_score, reduction = pre_embedding,
                                   do.fast = T, check_duplicates = check_duplicates)
  tSNE_distances_permuted = distances(pbmc.permuted@reductions$tsne@cell.embeddings)
  results.PCA <- list(reduced_dim_distances = tSNE_distances,
                      reduced_dim_distances_permuted = tSNE_distances_permuted)
  return(results.PCA)
}

Distances.UMAP_add_seu_obj <- function (pbmc, pbmc.permuted, K, pre_embedding = "pca", n = 30,
                                        m = 0.3, rerun = T)
{
  distances <- distances::distances
  if (rerun) {
    pbmc <- Seurat::RunUMAP(pbmc, dims = 1:K, seed.use = 100,
                            reduction = pre_embedding, n.neighbors = n, min.dist = m)
  }

  ## To obtain UMAP embeddings
  UMAP_data <- pbmc@reductions$umap@cell.embeddings
  ## write to the folder
  write_rds(UMAP_data, here::here(paste0("data/pbmc3k/pbmc_scdeed_umap_n_neighbors_", n, "_min_dist_", m, ".rds")))


  UMAP_distances = distances(UMAP_data)
  pbmc.permuted <- Seurat::RunUMAP(pbmc.permuted, dims = 1:K,
                                   seed.use = 100, reduction = pre_embedding, n.neighbors = n,
                                   min.dist = m)
  UMAP_distances_permuted = distances(pbmc.permuted@reductions$umap@cell.embeddings)
  results.PCA <- list(reduced_dim_distances = UMAP_distances,
                      reduced_dim_distances_permuted = UMAP_distances_permuted)
  return(results.PCA)
}

optimize_add_seu_obj <- function (input_data, input_data.permuted, pre_embedding, reduction.method,
                                  K, n, m, perplexity, results.PCA, similarity_percent, dubious_cutoff,
                                  trustworthy_cutoff, check_duplicates = T, rerun = T)
{
  if (reduction.method == "umap") {
    results <- Distances.UMAP(pbmc = input_data, pbmc.permuted = input_data.permuted,
                              K = K, pre_embedding = pre_embedding, n = n, m = m,
                              rerun = rerun)
  }
  else if (reduction.method == "tsne") {
    results <- Distances.tSNE_add_seu_obj(pbmc = input_data, pbmc.permuted = input_data.permuted,
                              K = K, pre_embedding = pre_embedding, perplexity_score = perplexity,
                              check_duplicates = check_duplicates, rerun = rerun)
  }
  similarity_score <- Cell.Similarity(results.PCA$pre_embedding_distances,
                                      results.PCA$pre_embedding_distances_permuted, results$reduced_dim_distances,
                                      results$reduced_dim_distances_permuted, similarity_percent)
  ClassifiedCells <- Cell.Classify(similarity_score$rho_original,
                                   similarity_score$rho_permuted, dubious_cutoff = dubious_cutoff,
                                   trustworthy_cutoff = trustworthy_cutoff)
  dub = ifelse(length(ClassifiedCells$dubious_cells) != 0,
               paste(ClassifiedCells$dubious_cells, sep = ",", collapse = ","),
               "none")
  int = ifelse(length(ClassifiedCells$intermediate_cells) !=
                 0, paste(ClassifiedCells$intermediate_cells, sep = ",",
                          collapse = ","), "none")
  trust = ifelse(length(ClassifiedCells$trustworthy_cells) !=
                   0, paste(ClassifiedCells$trustworthy_cells, sep = ",",
                            collapse = ","), "none")
  return(c(length(ClassifiedCells$dubious_cells), dub, trust,
           int))
}

scDEED_add_seu_obj <- function (input_data, K, n_neighbors = c(5, 20, 30, 40, 50),
                                min.dist = c(0.1, 0.4), similarity_percent = 0.5, reduction.method,
                                perplexity = c(seq(from = 20, to = 410, by = 30), seq(from = 450,
                                                                                      to = 800, by = 50)), pre_embedding = "pca", slot = "scale.data",
                                dubious_cutoff = 0.05, trustworthy_cutoff = 0.95, permuted = NA,
                                check_duplicates = T, rerun = T, default_assay = "active.assay")
{
  if (is.na(permuted)) {
    print("Permuting data")
    input_data.permuted <- suppressMessages(Permuted(input_data,
                                                     K = K, slot = slot,
                                                     default_assay = default_assay))
    print("Permutation finished")
  }
  else {
    input_data.permuted = permuted
  }
  results.PCA = suppressMessages(Distances.pre_embedding(input_data,
                                                         input_data.permuted,
                                                         K = K,
                                                         pre_embedding = pre_embedding))
  if (reduction.method == "umap") {
    all_pairs <- expand.grid(n_neighbors, min.dist)

    if (nrow(all_pairs) > 1) {
      print("Estimating time for each hyperparameter setting...")
      start = Sys.time()
      original = suppressMessages(foreach::`%do%`(foreach::foreach(n = all_pairs$Var1[1],
                                                                   m = all_pairs$Var2[1], .combine = "cbind"),
                                                  optimize(input_data, input_data.permuted, pre_embedding,
                                                           reduction.method, K, n = n, m = m, perplexity = NA,
                                                           results.PCA = results.PCA, similarity_percent = similarity_percent,
                                                           dubious_cutoff = dubious_cutoff, trustworthy_cutoff = trustworthy_cutoff,
                                                           rerun = rerun)))
      end = Sys.time()
      time = end - start
      print("Estimated time per hyperparameter setting ")
      print(time)
      print("Estimated time of completion")
      print(Sys.time() + time * nrow(all_pairs))
      original = as.matrix(original)
      all_dub <- suppressMessages(foreach::`%do%`(foreach::foreach(n = all_pairs$Var1[-1],
                                                                   m = all_pairs$Var2[-1], .combine = "cbind"),
                                                  optimize(input_data, input_data.permuted, pre_embedding,
                                                           reduction.method, K, n = n, m = m, perplexity = NA,
                                                           results.PCA = results.PCA, similarity_percent = similarity_percent,
                                                           dubious_cutoff = dubious_cutoff, trustworthy_cutoff = trustworthy_cutoff,
                                                           rerun = rerun)))
      if (nrow(all_pairs) == 2) {
        all_dub = as.matrix(all_dub)
      }
      all_dub = cbind(original, all_dub)
    }
    else {
      all_dub = suppressMessages(foreach::`%do%`(foreach::foreach(n = all_pairs$Var1,
                                                                  m = all_pairs$Var2, .combine = "cbind"), optimize(input_data,
                                                                                                                    input_data.permuted, pre_embedding, reduction.method,
                                                                                                                    K, n = n, m = m, perplexity = NA, results.PCA = results.PCA,
                                                                                                                    similarity_percent = similarity_percent, dubious_cutoff = dubious_cutoff,
                                                                                                                    trustworthy_cutoff = trustworthy_cutoff, rerun = rerun)))
      all_dub = as.matrix(all_dub)
    }
    all_dub = t(all_dub)
    colnames(all_dub) = c("number_dubious_cells", "dubious_cells",
                          "trustworthy_cells", "intermediate_cells")
    colnames(all_pairs) = c("n_neighbors", "min.dist")
    dubious_number_UMAP <- cbind(all_pairs, all_dub)
    dub_para <- data.frame(dubious_number_UMAP[, 1:3])
    dub_para_full = as.data.frame(dubious_number_UMAP)
  }
  if (reduction.method == "tsne") {
    if (is.null(input_data@reductions$tsne)) {
      input_data <- Seurat::RunTSNE(input_data, do.fast = T)
    }
    if (length(perplexity) > 1) {
      print("Estimating time for each hyperparameter setting...")
      start = Sys.time()
      original <- suppressMessages(foreach::`%do%`(foreach::foreach(p = perplexity[1],
                                                                    .combine = "cbind"), optimize_add_seu_obj(input_data, input_data.permuted,
                                                                                                  pre_embedding, reduction.method, K, n, m, perplexity = p,
                                                                                                  results.PCA = results.PCA, similarity_percent = similarity_percent,
                                                                                                  dubious_cutoff = dubious_cutoff, trustworthy_cutoff,
                                                                                                  check_duplicates = check_duplicates, rerun = rerun)))
      end = Sys.time()
      time = end - start
      print("Estimated time per hyperparameter setting ")
      print(time)
      original = as.matrix(original)
      print("Estimated time of completion")
      print(Sys.time() + time * length(perplexity))
      dubious_number_tSNE <- suppressMessages(foreach::`%do%`(foreach::foreach(p = perplexity[-1],
                                                                               .combine = "cbind"), optimize_add_seu_obj(input_data, input_data.permuted,
                                                                                                             pre_embedding, reduction.method, K, n, m, perplexity = p,
                                                                                                             results.PCA = results.PCA, similarity_percent = similarity_percent,
                                                                                                             dubious_cutoff = dubious_cutoff, trustworthy_cutoff,
                                                                                                             check_duplicates = check_duplicates, rerun = rerun)))
      if (length(perplexity) == 2) {
        dubious_number_tSNE = as.matrix(dubious_number_tSNE)
      }
      dubious_number_tSNE = cbind(original, dubious_number_tSNE)
    }
    else {
      dubious_number_tSNE <- suppressMessages(foreach::`%do%`(foreach::foreach(p = perplexity,
                                                                               .combine = "cbind"), optimize(input_data, input_data.permuted,
                                                                                                             pre_embedding, reduction.method, K, n, m, perplexity = p,
                                                                                                             results.PCA = results.PCA, similarity_percent = similarity_percent,
                                                                                                             dubious_cutoff = dubious_cutoff, trustworthy_cutoff,
                                                                                                             check_duplicates = check_duplicates, rerun = rerun)))
      dubious_number_tSNE = as.matrix(dubious_number_tSNE)
    }
    all_dub = t(dubious_number_tSNE)
    colnames(all_dub) = c("number_dubious_cells", "dubious_cells",
                          "trustworthy_cells", "intermediate_cells")
    dubious_number_tsne <- cbind(perplexity, all_dub)
    colnames(dubious_number_tsne)[1] = "perplexity"
    dub_para <- data.frame(dubious_number_tsne[, 1:2])
    if (length(dub_para) == 1) {
      dub_para = as.data.frame(t(dub_para))
    }
    dub_para_full = as.data.frame(dubious_number_tsne)
  }
  dub_para$number_dubious_cells = as.numeric(dub_para$number_dubious_cells)
  dub_para_full$number_dubious_cells = as.numeric(dub_para_full$number_dubious_cells)
  rownames(dub_para) = NULL
  rownames(dub_para_full) = NULL
  return(list(num_dubious = dub_para, full_results = dub_para_full))
}

#################################
load("data/pbmc3k/pbmc3k_clustering_UMAP.RData")
pbmc_umap <- umap.list$Festem@cell.embeddings |>
  tibble::as_tibble() |>
  mutate(cell_label = as.character(label.list$Festem)) |>
  rename(c("UMAP1" = "UMAP_1", "UMAP2" = "UMAP_2"))

#write_rds(pbmc_umap, "data/pbmc/pbmc_3k_festem/pbmc_umap.rds")

rownames(umap.list$Festem@cell.embeddings)

## Obtain gene names
rw_names <- rownames(umap.list$Festem@cell.embeddings) |> str_split("-")

name_vec <- c()
for (i in 1:length(rw_names)) {
  name <- rw_names[[i]][1]
  name_vec <- append(name_vec, name)

}

####################PCA
# library(RColorBrewer)
# my.color <- hue_pal()(13)
# names(my.color) <- 1:13
# umap.for.plot <- function(umap,cluster){
#   umap <- umap@cell.embeddings
#   umap <- as.data.frame(umap)
#   cbind(umap,cluster = cluster)
# }
load("data/pbmc3k/pbmc3k_hvggenes.RData")
pbmc <- readRDS("data/pbmc3k/pbmc3k_final.rds")
pbmc <- UpdateSeuratObject(pbmc)
# ref <- pbmc@active.ident
# names_fac <- names(pbmc@active.ident)
#
# ref <- as.character(ref)
# ref[c("CGGGCATGACCCAA","CTTGATTGATCTTC")] <- "Platelet"
# ref <- as.factor(ref)
# names(ref) <- names_fac

#pbmc <- pbmc[,!colnames(pbmc) %in% c("CGGGCATGACCCAA","CTTGATTGATCTTC")]
pbmc <- pbmc[,colnames(pbmc) %in% name_vec]

#pbmc <- pbmc[,ref!="Platelet"] ## Didn't run
gene.list <- list(EM[1:1000])
# ,
# hvgvst[1:1000],
# hvgdisp[1:1000],
# dub,
# devianceFS[1:1000],
# trendvar[1:1000]

# umap.list <- vector("list",length(gene.list))
# label.list <- vector("list",length(gene.list))
# names(label.list) <- c("Festem")
# names(umap.list) <- c("Festem")
# plots.list <- vector("list",length(gene.list))
for (i in 1:length(gene.list)){
  pbmc <- RunPCA(pbmc, verbose = FALSE,features = gene.list[[i]])
}

# valid_features <- gene.list[gene.list %in% rownames(pbmc)]
#
# pbmc <- subset(pbmc, features = valid_features)

## For tSNE
result <- scDEED_add_seu_obj(pbmc, K = 9, reduction.method = 'tsne', perplexity = c(5, 15, 30, 51))
result$num_dubious
saveRDS(result$num_dubious, 'data/pbmc3k/pbmc_scdeed_tsne_results.rds')

## For UMAP
result_umap <- scDEED_add_seu_obj(pbmc, K = 9,
                     reduction.method = 'umap',
                     n_neighbors = c(5, 15, 30),
                     min.dist = c(0.01, 0.3, 0.99))

#head(result$num_dubious)
result_umap$num_dubious
saveRDS(result_umap$num_dubious, 'data/pbmc3k/pbmc_scdeed_umap_results.rds')
