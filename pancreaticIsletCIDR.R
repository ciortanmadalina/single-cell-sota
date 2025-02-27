## install CIDR from https://github.com/VCCRI/CIDR
## install the package plot3D for 3D plots
## install the package rgl for interactive 3D plots

"
Human Pancreatic Islet scRNA-Seq Dataset
CIDR-examples contains a human pancreatic islet single-cell RNA-Seq dataset, located in the PancreaticIslet folder. 
In this dataset there are 60 cells in 6 cell types after we exclude undefined cells and bulk RNA-Seq samples.

Reference for the human pancreatic islet dataset:

Li, J. et al. Single-cell transcriptomes reveal characteristic features of human pancreatic islet cell types. 
EMBO Reports 17, 178–187 (2016).
"
rm(list = ls())
.rs.restartR()

library(cidr)

## Read in data 
## Tag tables were downloaded from the data repository NCBI Gene Expression Omnibus (GSE73727)
## Tag tables were combined into one table and undefined cells and bulk RNA-Seq samples have been excluded.
## Pseudo-tag 1 has been subtracted from each entry.
pancreaticIsletTags <- read.csv("single-cell-sota/input/pancreaticIsletCIDR/pancreaticIsletTags.csv")
rownames(pancreaticIsletTags) <- pancreaticIsletTags[,1]
pancreaticIsletTags <- pancreaticIsletTags[,-1]
info <- read.csv("single-cell-sota/input/pancreaticIsletCIDR/SraRunTable.txt",sep="\t")
cellType <- info$assigned_cell_type_s[match(colnames(pancreaticIsletTags),info$Sample_Name_s)]

truth <- info$assigned_cell_type_s[match(colnames(pancreaticIsletTags),info$Sample_Name_s)]
truth <- data.frame(colnames(pancreaticIsletTags) , truth)
write.csv(truth,'single-cell-sota/input/pancreaticIsletCIDR/truth.csv')

cellType <- factor(cellType)
types <- levels(cellType)

## Assign each cell type a color
scols <- c("red","green","blue","grey","purple","brown")
cols <- rep(NA,length(cellType))
for (i in 1:length(cols)){
  cols[i] <- scols[which(types==cellType[i])]
}

# Standard principal component analysis using prcomp
priorTPM <- 1
pan10 <- pancreaticIsletTags[rowSums(pancreaticIsletTags)>10,]
pan10_lcpm <- log2(t(t(pan10)/colSums(pan10))*1000000+priorTPM)
PC_lcpm <- prcomp(t(pan10_lcpm))
plot(PC_lcpm$x[,c(1,2)], col=cols, pch=1, 
     xlab="PC1", ylab="PC2", 
     main="Principal Component Analysis (prcomp)")
legend("topright", legend=types, col=scols, pch=1)

###########
## CIDR ###
###########
scPan <- scDataConstructor(as.matrix(pancreaticIsletTags))
scPan <- determineDropoutCandidates(scPan)
scPan <- wThreshold(scPan)
scPan <- scDissim(scPan)
scPan <- scPCA(scPan)
scPan <- nPC(scPan)
scPan <- scCluster(scPan)

## Two dimensional visualization plot output by CIDR
## Different colors denote the cell types annotated by the human pancreatic islet single-cell RNA-Seq study
## Different plotting symbols denote the clusters output by CIDR
plot(scPan@PC[,c(1,2)], col=cols, pch=scPan@clusters, 
     main="CIDR", xlab="PC1", ylab="PC2")
legend("bottomleft", legend=types, col=scols, pch=15)

## 3D
library(plot3D)
scatter3D(scPan@PC[,1], scPan@PC[,2], scPan@PC[,3],
          xlab="PC1", ylab="PC2", zlab="PC3",
          colvar=NULL, col=cols, pch=scPan@clusters, 
          phi=5, theta=55, pch=20)
legend("bottomright", legend=types, col=scols, pch=15)

## 3D - interactive
library(rgl)
plot3d(scPan@PC[,1:3], col=cols, 
       xlab="PC1", ylab="PC2", zlab="PC3")

## Use Adjusted Rand Index to measure the accuracy of CIDR clustering
ARI_CIDR <- adjustedRandIndex(scPan@clusters,cols)
ARI_CIDR
##  0.6830087

