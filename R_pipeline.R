rm(list = ls())
.rs.restartR()

library(Matrix)
library(SingleCellExperiment)
library(scater)
library(pcaMethods)
library(pcaReduce)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
options(stringsAsFactors = FALSE)

df<-get(load('all_c_elegans/data/NeuronalGeneCount'))
rm(NeuronalGenCount)
df <- as.matrix(df)
anno <- colnames(df)
head(df[ , 1:3])

umi <- SingleCellExperiment(
  assays = list(counts = as.matrix(df)), 
  colData = anno
)
dim(umi)

checkPrint <- function(x, cellidx, isCounts ){
  if (isCounts == FALSE){
    o = which(counts(x)!=0,arr.ind = T)
    #cellidx = 2
    counts(x)[head(unname(o[o[, "col"] == cellidx, 'row']), 15),cellidx]
  }
  else {
    o = which(x!=0,arr.ind = T)
    #cellidx = 2
    x[head(unname(o[o[, "col"] == cellidx, 'row']), 15),cellidx]
  }
}
checkPrint(umi, 2, FALSE)

keep_feature <- rowSums(counts(umi) > 0) > 0
umi <- umi[keep_feature, ]

# define feature names in feature_symbol column
rowData(umi)$feature_symbol <- rownames(umi)
# remove features with duplicated names
umi <- umi[!duplicated(rowData(umi)$feature_symbol), ]

umi <- calculateQCMetrics(
  umi
)

hist(
  umi$total_counts,
  breaks = 100
)
abline(v = 10, col = "red")

filter_by_total_counts <- (umi$total_counts > 100)
table(filter_by_total_counts)

hist(
  umi$total_features,
  breaks = 100
)
abline(v = 50, col = "red")


filter_by_expr_features <- (umi$total_features > 50)
table(filter_by_expr_features)

umi$use <- (
  # sufficient features (genes)
  filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts 
)
table(umi$use)
#Automatic filtering
umi <- plotPCA(
  umi,
  size_by = "total_features", 
  shape_by = "use",
  pca_data_input = "pdata",
  detect_outliers = TRUE,
  return_SCE = TRUE
)

table(umi$outlier)

library(limma)
auto <- colnames(umi)[umi$outlier]
man <- colnames(umi)[!umi$use]
venn.diag <- vennCounts(
  cbind(colnames(umi) %in% auto,
        colnames(umi) %in% man)
)
vennDiagram(
  venn.diag,
  names = c("Automatic", "Manual"),
  circle.col = c("blue", "green")
)

# combine outliers
umi$use <- (
  umi$use  &
    !umi$outlier
)


filter_genes <- apply(
  counts(umi[ , colData(umi)$use]), 
  1, 
  function(x) length(x[x > 1]) >= 2
)
rowData(umi)$use <- filter_genes

table(umi$use)

assay(umi, "logcounts_raw") <- log2(counts(umi) + 1)
reducedDim(umi) <- NULL

umi.qc <- umi[rowData(umi)$use, colData(umi)$use]

umi$qc
head(umi$rowData[ , 1:3])
umi

calc_cpm <-  function (expr_mat, spikes = NULL) 
  {
    norm_factor <- colSums(expr_mat[-spikes, ])
    return(t(t(expr_mat)/norm_factor)) * 10^6
}

checkPrint(umi, 2, FALSE)
logcounts(umi.qc) <- log2(calculateCPM(umi.qc, use.size.factors = FALSE) + 1)
checkPrint(logcounts(umi.qc), 2, TRUE)
checkPrint(logcounts_raw(umi.qc), 2, TRUE)

# https://bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/SC3.html
umi.qc <- sc3_estimate_k(umi.qc)
metadata(umi.qc)$sc3$k_estimation
plotPCA(umi.qc)

umi
umi.qc <- sc3(umi.qc, ks = 50, biology = TRUE)

col_data <- colData(umi)
head(col_data[ , grep("sc3_", colnames(col_data))])

plotPCA(
  umi, 
  colour_by = "sc3_3_clusters", 
  size_by = "sc3_3_log2_outlier_score"
)

sc3_plot_consensus(umi, k = 118, show_pdata = "cell_type2")

sc3_plot_silhouette(umi, k = 118)



# PCA reduce
input <- logcounts(umi.qc[rowData(umi.qc)$sc3_gene_filter, ])
# run pcaReduce 1 time creating hierarchies from 1 to 30 clusters
pca.red <- PCAreduce(t(input), nbt = 1, q = 30, method = 'S')[[1]]

colData(umi.qc)$pcaReduce <- as.character(pca.red[,32 - 10])
plotPCA(umi.qc, colour_by = "pcaReduce")

colData(umi.qc)$pcaReduce <- as.character(pca.red[,32 - 2])
plotPCA(umi.qc, colour_by = "pcaReduce")

colData(umi.qc)$pcaReduce <- as.character(pca.red[,32 - 10])
adjustedRandIndex(colData(umi.qc)$cell_type2, colData(deng)$pcaReduce)

## tSNE + kmeans
umi <- plotTSNE(umi.qc, rand_seed = 1, return_SCE = TRUE)

colData(umi)$tSNE_kmeans <- as.character(kmeans(umi.@reducedDims$TSNE, centers = 8)$clust)
plotTSNE(umi.qc, rand_seed = 1, colour_by = "tSNE_kmeans")

colData(umi)$tSNE_kmeans <- as.character(kmeans(umi@reducedDims$TSNE, centers = 10)$clust)
#adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$tSNE_kmeans)

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

colData(umi)$SNNCliq <- as.character(snn.res[,1])
plotPCA(umi, colour_by = "SNNCliq")

#adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$SNNCliq)

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
