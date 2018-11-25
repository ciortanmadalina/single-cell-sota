rm(list = ls())
.rs.restartR()
library(pcaMethods)
library(pcaReduce)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
set.seed(1234567)
deng <- readRDS("deng/deng-reads.rds")

## Write as csv
write.csv(counts(deng),'single-cell-sota/input/deng/deng.csv')
write.csv(colData(deng)$cell_type2,'single-cell-sota/input/deng/truth.csv')


dim(deng)
table(colData(deng)$cell_type2)

plotPCA(deng, colour_by = "cell_type2")
deng <- sc3_estimate_k(deng)
metadata(deng)$sc3$k_estimation
plotPCA(deng, colour_by = "cell_type1")
deng <- sc3(deng, ks = 8:10, biology = TRUE)
sc3_plot_consensus(deng, k = 10, show_pdata = "cell_type2")
sc3_plot_silhouette(deng, k = 10)
adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$sc3_10_clusters)

sc3_interactive(deng)

# PCA reduce
input <- logcounts(deng[rowData(deng)$sc3_gene_filter, ])
# run pcaReduce 1 time creating hierarchies from 1 to 30 clusters
pca.red <- PCAreduce(t(input), nbt = 1, q = 30, method = 'S')[[1]]

colData(deng)$pcaReduce <- as.character(pca.red[,32 - 10])
plotPCA(deng, colour_by = "pcaReduce")

colData(deng)$pcaReduce <- as.character(pca.red[,32 - 2])
plotPCA(deng, colour_by = "pcaReduce")

colData(deng)$pcaReduce <- as.character(pca.red[,32 - 10])
adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$pcaReduce)

## tSNE + kmeans
deng <- plotTSNE(deng, rand_seed = 1, return_SCE = TRUE)

colData(deng)$tSNE_kmeans <- as.character(kmeans(deng@reducedDims$TSNE, centers = 8)$clust)
plotTSNE(deng, rand_seed = 1, colour_by = "tSNE_kmeans")

colData(deng)$tSNE_kmeans <- as.character(kmeans(deng@reducedDims$TSNE, centers = 10)$clust)
adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$tSNE_kmeans)

##SNN-Cliq
distan <- "euclidean"
par.k <- 3
par.r <- 0.7
par.m <- 0.5
# construct a graph
scRNA.seq.funcs::SNN(
  data = t(input),
  outfile = "snn-cliq.txt",
  k = par.k,
  distance = distan
)
# find clusters in the graph
snn.res <- 
  system(
    paste0(
      "python utils/Cliq.py ", 
      "-i snn-cliq.txt ",
      "-o res-snn-cliq.txt ",
      "-r ", par.r,
      " -m ", par.m
    ),
    intern = TRUE
  )
cat(paste(snn.res, collapse = "\n"))
snn.res <- read.table("res-snn-cliq.txt")
# remove files that were created during the analysis
system("rm snn-cliq.txt res-snn-cliq.txt")

colData(deng)$SNNCliq <- as.character(snn.res[,1])
plotPCA(deng, colour_by = "SNNCliq")

adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$SNNCliq)

## SINCERA
# perform gene-by-gene per-sample z-score transformation
dat <- apply(input, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
# hierarchical clustering
dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
hc <- hclust(dd, method = "average")


num.singleton <- 0
kk <- 1
for (i in 2:dim(dat)[2]) {
  clusters <- cutree(hc, k = i)
  clustersizes <- as.data.frame(table(clusters))
  singleton.clusters <- which(clustersizes$Freq < 2)
  if (length(singleton.clusters) <= num.singleton) {
    kk <- i
  } else {
    break;
  }
}
cat(kk)

pheatmap(
  t(dat),
  cluster_cols = hc,
  cutree_cols = kk,
  kmeans_k = 100,
  show_rownames = FALSE
)

colData(deng)$SINCERA <- as.character(cutree(hc, k = kk))
adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$SINCERA)


## Feature selection 

library(scRNA.seq.funcs)
library(matrixStats)
library(M3Drop)
library(RColorBrewer)
library(SingleCellExperiment)
set.seed(1)

deng <- readRDS("deng/deng-reads.rds")
cellLabels <- colData(deng)$cell_type2


deng_list <- M3DropCleanData(
  counts(deng),
  labels = cellLabels,
  min_detected_genes = 100,
  is.counts = TRUE
)
expr_matrix <- deng_list$data # Normalized & filtered expression matrix
celltype_labs <- factor(deng_list$labels) # filtered cell-type labels
cell_colors <- brewer.pal(max(3,length(unique(celltype_labs))), "Set3")



#How many cells & genes have been removed by this filtering? 
ncol(counts(deng)) - ncol(deng_list$data)


Brennecke_HVG <- BrenneckeGetVariableGenes(
  expr_matrix,
  fdr = 0.01,
  minBiolDisp = 0.5
)
HVG_genes <- Brennecke_HVG$Gene

# High dropout genes

K <- 49
S_sim <- 10^seq(from = -3, to = 4, by = 0.05) # range of expression values
MM <- 1 - S_sim / (K + S_sim)
plot(
  S_sim, 
  MM, 
  type = "l", 
  lwd = 3, 
  xlab = "Expression", 
  ylab = "Dropout Rate", 
  xlim = c(1,1000)
)
S1 <- 10
P1 <- 1 - S1 / (K + S1) # Expression & dropouts for cells in condition 1
S2 <- 750
P2 <- 1 - S2 / (K + S2) # Expression & dropouts for cells in condition 2
points(
  c(S1, S2),
  c(P1, P2), 
  pch = 16, 
  col = "grey85", 
  cex = 3
)
mix <- 0.5 # proportion of cells in condition 1
points(
  S1 * mix + S2 * (1 - mix), 
  P1 * mix + P2 * (1 - mix), 
  pch = 16, 
  col = "grey35", 
  cex = 3
)


M3Drop_genes <- M3DropFeatureSelection(
  expr_matrix,
  mt_method = "fdr",
  mt_threshold = 0.01
)

M3Drop_genes <- M3Drop_genes$Gene


deng_int <- NBumiConvertToInteger(counts(deng))
DANB_fit <- NBumiFitModel(deng_int) # DANB is fit to the raw count matrix
# Perform DANB feature selection
DropFS <- NBumiFeatureSelectionCombinedDrop(DANB_fit)
DANB_genes <- names(DropFS[1:1500])

#Correlated expression
cor_mat <- cor(t(expr_matrix), method = "spearman") # Gene-gene correlations
diag(cor_mat) <- rep(0, times = nrow(expr_matrix))
score <- apply(cor_mat, 1, function(x) {max(abs(x))}) #Correlation of highest magnitude
names(score) <- rownames(expr_matrix);
score <- score[order(-score)]
Cor_genes <- names(score[1:1500])

#PCA
# PCA is typically performed on log-transformed expression data
pca <- prcomp(log(expr_matrix + 1) / log(2))

# plot projection
plot(
  pca$rotation[,1], 
  pca$rotation[,2], 
  pch = 16, 
  col = cell_colors[as.factor(celltype_labs)]
) 

score <- rowSums(abs(pca$x[,c(1,2)])) 
names(score) <- rownames(expr_matrix)
score <- score[order(-score)]
PCA_genes <- names(score[1:1500])


# Compare methods
M3DropExpressionHeatmap(
  M3Drop_genes,
  expr_matrix,
  cell_labels = celltype_labs
)

J <- sum(M3Drop_genes %in% HVG_genes)/length(unique(c(M3Drop_genes, HVG_genes)))
