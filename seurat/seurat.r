rm(list = ls())
.rs.restartR()

library(SingleCellExperiment)
library(Seurat)
library(mclust)
library(dplyr)
options(stringsAsFactors = FALSE)
# Pipeline inspired by
# https://hemberg-lab.github.io/scRNA.seq.course/seurat-chapter.html

# Load deng annotated data
deng <- readRDS("single-cell-sota/input/deng-reads.rds")

# Or load c-elegans data
df<-get(load('single-cell-sota/input/NeuronalGeneCount'))
rm(NeuronalGenCount)
df <- as.matrix(df)
anno <- colnames(df)
head(df[ , 1:3])
deng <- SingleCellExperiment(
  assays = list(counts = as.matrix(df)), 
  colData = anno
)
rm(df) # remove heavy object

deng # check single cell experiment object

seuset <- CreateSeuratObject(
  raw.data = counts(deng),
  min.cells = 3, 
  min.genes = 200
)

VlnPlot(
  object = seuset, 
  features.plot = c("nGene", "nUMI"), 
  nCol = 2
)

GenePlot(
  object = seuset, 
  gene1 = "nUMI", 
  gene2 = "nGene"
)

seuset <- FilterCells(
  object = seuset, 
  subset.names = c("nUMI"), 
  high.thresholds = c(4000)#c(2e7)
)


seuset <- NormalizeData(
  object = seuset, 
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)

seuset <- FindVariableGenes(
  object = seuset,
  mean.function = ExpMean, 
  dispersion.function = LogVMR, 
  x.low.cutoff = 0.0125, 
  x.high.cutoff = 3, 
  y.cutoff = 0.5
)

length(x = seuset@var.genes)

# Confounders
seuset <- ScaleData(
  object = seuset, 
  vars.to.regress = c("nUMI")
)

seuset <- RunPCA(
  object = seuset, 
  pc.genes = seuset@var.genes, 
  do.print = TRUE, 
  pcs.print = 1:5, 
  genes.print = 5
)

PrintPCA(object = seuset, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

VizPCA(object = seuset, pcs.use = 1:2)

PCAPlot(object = seuset, dim.1 = 1, dim.2 = 2)

PCHeatmap(
  object = seuset, 
  pc.use = 1:6, 
  cells.use = 500, 
  do.balanced = TRUE, 
  label.columns = FALSE,
  use.full = FALSE
)

seuset <- JackStraw(
  object = seuset, 
  num.replicate = 100, 
  do.print = FALSE
)


JackStrawPlot(object = seuset, PCs = 1:9)

PCElbowPlot(object = seuset)

seuset <- FindClusters(
  object = seuset, 
  reduction.type = "pca", 
  dims.use = 1:8, 
  resolution = 1.0, 
  print.output = 0, 
  save.SNN = TRUE
)

PrintFindClustersParams(object = seuset)


table(seuset@ident)
write.csv(seuset@ident,'single-cell-sota/output/r_elegans_seurat.csv')


# Validate against ground truth
adjustedRandIndex(colData(deng)[seuset@cell.names, ]$cell_type2, seuset@ident)

seuset <- RunTSNE(
  object = seuset,
  dims.use = 1:8,
  do.fast = TRUE
)
TSNEPlot(object = seuset)

seuset <- RunTSNE(
  object = seuset,
  dims.use = 1:8,
  do.fast = TRUE
)
TSNEPlot(object = seuset)

markers2 <- FindMarkers(seuset, 2)
VlnPlot(object = seuset, features.plot = rownames(markers2)[1:2])

FeaturePlot(
  seuset, 
  head(rownames(markers2)), 
  cols.use = c("lightgrey", "blue"), 
  nCol = 3
)
