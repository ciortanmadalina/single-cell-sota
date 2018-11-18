#https://davetang.org/muse/2017/10/01/getting-started-monocle/
rm(list = ls())
.rs.restartR()

library(monocle)
#source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")
library(cellrangerRkit)
my_dir <- "single-cell-sota/input"
# load data
gbm <- load_cellranger_matrix(my_dir)
class(gbm)
dim(exprs(gbm))
exprs(gbm)[1:5, 1:5]
# check out the phenotypic data
dim(pData(gbm))
head(pData(gbm))

# check out the feature information
dim(fData(gbm))
# this table will be useful for matching gene IDs to symbols
head(fData(gbm))

# warning!
my_cds <- newCellDataSet(exprs(gbm),
                         phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                         featureData = new("AnnotatedDataFrame", data = fData(gbm)),
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
# rename gene symbol column
my_feat <- fData(gbm)
names(my_feat) <- c('id', 'gene_short_name')
# no warning
my_cds <- newCellDataSet(exprs(gbm),
                         phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                         featureData = new("AnnotatedDataFrame", data = my_feat),
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

my_cds

# normalisation and variance estimation steps, which will be used in the differential expression analyses later on.
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds) # freezes computer for 5 min

my_cds <- detectGenes(my_cds, min_expr = 0.1)
head(fData(my_cds))

sum((exprs(my_cds['ENSG00000239945',])))
sum((exprs(my_cds['ENSG00000238009',])))
head(pData(my_cds))
sum((exprs(my_cds)[,"AAACCTGAGCATCATC-1"])>0)
summary(pData(my_cds)$num_genes_expressed)

# standardise to Z-distribution
x <- pData(my_cds)$num_genes_expressed
x_1 <- (x - mean(x)) / sd(x)
summary(x_1)


library(ggplot2)
# I like the default theme of cowplot
library(cowplot)

df <- data.frame(x = x_1)
ggplot(df, aes(x)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')

# add a UMI column into phenoData
# use Matrix because data is in the sparse format
pData(my_cds)$UMI <- Matrix::colSums(exprs(my_cds))
head(pData(my_cds))

ggplot(pData(my_cds), aes(num_genes_expressed, UMI)) + geom_point()

# select gene based on their average expression and variability across cells. 
#The dispersionTable() function calculates the mean and dispersion values.
disp_table <- dispersionTable(my_cds)
head(disp_table)

table(disp_table$mean_expression>=0.1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
dim(unsup_clustering_genes)

my_cds <- setOrderingFilter(my_cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(my_cds)

plot_pc_variance_explained(my_cds, return_all = FALSE)

# if num_dim is not set, the default is set to 50
# I'll use 5 here based on plot_pc_variance_explained
my_cds <- reduceDimension(my_cds, max_components = 2, num_dim = 5,
                          reduction_method = 'tSNE', verbose = TRUE)


# perform unsupervised clustering requesting 15-1 clusters
my_cds <- clusterCells(my_cds, num_clusters = 15)

# cluster information is in the Cluster column
head(pData(my_cds))

# store cluster info for comparing later
my_cluster_dim_5 <- pData(my_cds)$Cluster

plot_cell_clusters(my_cds)

# Try 10 clusters 
my_cds <- reduceDimension(my_cds, max_components = 2, num_dim = 10,
reduction_method = 'tSNE', verbose = TRUE)
my_cds <- clusterCells(my_cds, num_clusters = 15)

my_cluster_dim_10 <- pData(my_cds)$Cluster
plot_cell_clusters(my_cds)

# Validate results with adjusted rand index
library(clues)
adjustedRand(as.numeric(my_cluster_dim_5), as.numeric(my_cluster_dim_10))

my_cds <- clusterCells(my_cds)
head(pData(my_cds))
