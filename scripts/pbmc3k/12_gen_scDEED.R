library(scDEED)
result = scDEED(pbmc, K = 9, reduction.method = 'tsne', perplexity = c(5, 15, 30, 40, 50))
result_umap = scDEED(pbmc, K = 9, reduction.method = 'umap')

head(result$num_dubious)
head(result_umap$num_dubious)
