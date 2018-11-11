rm(list = ls())
.rs.restartR()

# Cleaning expression matrix

#To illustrate cell QC, we consider a dataset of induced pluripotent stem cells generated from three 
#different individuals (Tung et al. 2017) in Yoav Gilad’s lab at the University of Chicago. 
#The experiments were carried out on the Fluidigm C1 platform and to facilitate the quantification both 
#unique molecular identifiers (UMIs) and ERCC spike-ins were used. The data files are located in the tung folder 
#in your working directory. These files are the copies of the original files made on the 15/03/16. 
#We will use these copies for reproducibility purposes.


library(Matrix)
library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)
# https://hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html#exprs-qc
molecules <- read.table("tung/molecules.txt", sep = "\t")
anno <- read.table("tung/annotation.txt", sep = "\t", header = TRUE)

umi <- SingleCellExperiment(
  assays = list(counts = as.matrix(molecules)), 
  colData = anno
)

# Remove genes that are not expressed in any cell
keep_feature <- rowSums(counts(umi) > 0) > 0
umi <- umi[keep_feature, ]

# Define control features (genes) - ERCC spike-ins and mitochondrial genes (provided by the authors):

isSpike(umi, "ERCC") <- grepl("^ERCC-", rownames(umi))
isSpike(umi, "MT") <- rownames(umi) %in% 
  c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
    "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
    "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
    "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
    "ENSG00000198840")
# Calculate the quality metrics: this adds total_counts, etc

umi <- calculateQCMetrics(
  umi,
  feature_controls = list(
    ERCC = isSpike(umi, "ERCC"), 
    MT = isSpike(umi, "MT")
  )
)

hist(
  umi$total_counts,
  breaks = 100
)
abline(v = 25000, col = "red")

filter_by_total_counts <- (umi$total_counts > 25000)
table(filter_by_total_counts)

hist(
  umi$total_features,
  breaks = 100
)
abline(v = 7000, col = "red")


filter_by_expr_features <- (umi$total_features > 7000)
table(filter_by_expr_features)

plotPhenoData(
  umi,
  aes_string(
    x = "total_features",
    y = "pct_counts_MT",
    colour = "batch"
  )
)

plotPhenoData(
  umi,
  aes_string(
    x = "total_features",
    y = "pct_counts_ERCC",
    colour = "batch"
  )
)

filter_by_ERCC <- umi$batch != "NA19098.r2"
table(filter_by_ERCC)
filter_by_MT <- umi$pct_counts_MT < 10
table(filter_by_MT)

# Cell filtering

umi$use <- (
  # sufficient features (genes)
  filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # sufficient endogenous RNA
    filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    filter_by_MT
)
table(umi$use)


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

# Gene analysis
plotQC(umi, type = "highest-expression")

filter_genes <- apply(
  counts(umi[ , colData(umi)$use]), 
  1, 
  function(x) length(x[x > 1]) >= 2
)
rowData(umi)$use <- filter_genes
table(filter_genes)

dim(umi[rowData(umi)$use, colData(umi)$use])

assay(umi, "logcounts_raw") <- log2(counts(umi) + 1)
reducedDim(umi) <- NULL
saveRDS(umi, file = "all_c_elegans/data/umi.rds")

library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)
umi <- readRDS("all_c_elegans/data/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control


plotPCA(
  umi[endog_genes, ],
  exprs_values = "counts",
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)

plotPCA(
  umi[endog_genes, ],
  exprs_values = "logcounts_raw",
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)

plotPCA(
  umi.qc[endog_genes, ],
  exprs_values = "logcounts_raw",
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)

plotTSNE(
  umi[endog_genes, ],
  exprs_values = "logcounts_raw",
  perplexity = 130,
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  rand_seed = 123456
)

plotTSNE(
  umi.qc[endog_genes, ],
  exprs_values = "logcounts_raw",
  perplexity = 130,
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  rand_seed = 123456
)
# Cofounding factors


library(scater, quietly = TRUE)
options(stringsAsFactors = FALSE)
umi <- readRDS("all_c_elegans/data/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control


plotPCA(
  umi.qc[endog_genes, ],
  exprs_values = "logcounts_raw",
  colour_by = "batch",
  size_by = "total_features"
)

plotQC(
  umi.qc[endog_genes, ],
  type = "find-pcs",
  exprs_values = "logcounts_raw",
  variable = "total_features"
)


plotQC(
  umi.qc[endog_genes, ],
  type = "expl",
  exprs_values = "logcounts_raw",
  variables = c(
    "total_features",
    "total_counts",
    "batch",
    "individual",
    "pct_counts_ERCC",
    "pct_counts_MT"
  )
)
# Normalization
# correct for library size by multiplying with total count of genes by cell

calc_cpm <-
  function (expr_mat, spikes = NULL) 
  {
    norm_factor <- colSums(expr_mat[-spikes, ])
    return(t(t(expr_mat)/norm_factor)) * 10^6
  }

# size factor First the geometric mean of each gene across all cells is calculated. The size factor 
#for each cell is the median across genes of the ratio of the expression to the gene’s geometric mean.
calc_sf <-
  function (expr_mat, spikes = NULL) 
  {
    geomeans <- exp(rowMeans(log(expr_mat[-spikes, ])))
    SF <- function(cnts) {
      median((cnts/geomeans)[(is.finite(geomeans) & geomeans > 
                                0)])
    }
    norm_factor <- apply(expr_mat[-spikes, ], 2, SF)
    return(t(t(expr_mat)/norm_factor))
  }

#The upperquartile (UQ) was proposed by (Bullard et al. 2010). Here each column is divided by the 75% quantile
#of the counts for each library. Often the calculated quantile is scaled by the median across cells to keep 
#the absolute level of expression relatively consistent. 


calc_uq <-
  function (expr_mat, spikes = NULL) 
  {
    UQ <- function(x) {
      quantile(x[x > 0], 0.75)
    }
    uq <- unlist(apply(expr_mat[-spikes, ], 2, UQ))
    norm_factor <- uq/median(uq)
    return(t(t(expr_mat)/norm_factor))
  }


Down_Sample_Matrix <-
  function (expr_mat) 
  {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }


calc_cell_RLE <-
  function (expr_mat, spikes = NULL) 
  {
    RLE_gene <- function(x) {
      if (median(unlist(x)) > 0) {
        log((x + 1)/(median(unlist(x)) + 1))/log(2)
      }
      else {
        rep(NA, times = length(x))
      }
    }
    if (!is.null(spikes)) {
      RLE_matrix <- t(apply(expr_mat[-spikes, ], 1, RLE_gene))
    }
    else {
      RLE_matrix <- t(apply(expr_mat, 1, RLE_gene))
    }
    cell_RLE <- apply(RLE_matrix, 2, median, na.rm = T)
    return(cell_RLE)
  }


library(scRNA.seq.funcs)
library(scater)
library(scran)
options(stringsAsFactors = FALSE)
set.seed(1234567)
umi <- readRDS("all_c_elegans/data/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control

plotPCA(
  umi.qc[endog_genes, ],
  exprs_values = "logcounts_raw",
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)


logcounts(umi.qc) <- log2(calculateCPM(umi.qc, use.size.factors = FALSE) + 1)
plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)
plotRLE(
  umi.qc[endog_genes, ], 
  exprs_mats = list(Raw = "logcounts_raw", CPM = "logcounts"),
  exprs_logged = c(TRUE, TRUE),
  colour_by = "batch"
)


umi.qc <- normaliseExprs(
  umi.qc,
  method = "RLE", 
  feature_set = endog_genes,
  return_log = TRUE,
  return_norm_as_exprs = TRUE
)
## Warning in normalizeSCE(object, exprs_values = exprs_values, return_log
## = return_log, : spike-in transcripts in 'ERCC' should have their own size
## factors
plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)

plotRLE(
  umi.qc[endog_genes, ], 
  exprs_mats = list(Raw = "logcounts_raw", RLE = "logcounts"),
  exprs_logged = c(TRUE, TRUE),
  colour_by = "batch"
)


umi.qc <- normaliseExprs(
  umi.qc,
  method = "upperquartile", 
  feature_set = endog_genes,
  p = 0.99,
  return_log = TRUE,
  return_norm_as_exprs = TRUE
)
## Warning in normalizeSCE(object, exprs_values = exprs_values, return_log
## = return_log, : spike-in transcripts in 'ERCC' should have their own size
## factors
plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)

plotRLE(
  umi.qc[endog_genes, ], 
  exprs_mats = list(Raw = "logcounts_raw", UQ = "logcounts"),
  exprs_logged = c(TRUE, TRUE),
  colour_by = "batch"
)


umi.qc <- normaliseExprs(
  umi.qc,
  method = "TMM",
  feature_set = endog_genes,
  return_log = TRUE,
  return_norm_as_exprs = TRUE
)
## Warning in normalizeSCE(object, exprs_values = exprs_values, return_log
## = return_log, : spike-in transcripts in 'ERCC' should have their own size
## factors
plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)

plotRLE(
  umi.qc[endog_genes, ], 
  exprs_mats = list(Raw = "logcounts_raw", TMM = "logcounts"),
  exprs_logged = c(TRUE, TRUE),
  colour_by = "batch"
)


qclust <- quickCluster(umi.qc, min.size = 30)
umi.qc <- computeSumFactors(umi.qc, sizes = 15, clusters = qclust)
umi.qc <- normalize(umi.qc)
## Warning in .local(object, ...): spike-in transcripts in 'ERCC' should have
## their own size factors
plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)

plotRLE(
  umi.qc[endog_genes, ], 
  exprs_mats = list(Raw = "logcounts_raw", scran = "logcounts"),
  exprs_logged = c(TRUE, TRUE),
  colour_by = "batch"
)

umi.qc <- getBMFeatureAnnos(
  umi.qc,
  filters = "ensembl_gene_id", 
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol",
    "chromosome_name",
    "start_position",
    "end_position"
  ), 
  feature_symbol = "hgnc_symbol",
  feature_id = "ensembl_gene_id",
  biomart = "ENSEMBL_MART_ENSEMBL", 
  dataset = "hsapiens_gene_ensembl",
  host = "www.ensembl.org"
)

umi.qc.ann <- umi.qc[!is.na(rowData(umi.qc)$ensembl_gene_id), ]

eff_length <- 
  abs(rowData(umi.qc.ann)$end_position - rowData(umi.qc.ann)$start_position) / 1000

plot(eff_length, rowMeans(counts(umi.qc.ann)))

tpm(umi.qc.ann) <- log2(calculateTPM(umi.qc.ann, eff_length) + 1)

plotPCA(
  umi.qc.ann,
  exprs_values = "tpm",
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)
tpm(umi.qc.ann) <- log2(calculateFPKM(umi.qc.ann, eff_length) + 1)
plotPCA(
  umi.qc.ann,
  exprs_values = "tpm",
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)

## COFOUNDERS

library(scRNA.seq.funcs)
library(RUVSeq)
library(scater)
library(SingleCellExperiment)
library(scran)
library(kBET)
library(sva) # Combat
library(edgeR)
set.seed(1234567)
options(stringsAsFactors = FALSE)
umi <- readRDS("all_c_elegans/data/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control
erccs <- rowData(umi.qc)$is_feature_control

qclust <- quickCluster(umi.qc, min.size = 30)
umi.qc <- computeSumFactors(umi.qc, sizes = 15, clusters = qclust)
umi.qc <- normalize(umi.qc)


ruvg <- RUVg(counts(umi.qc), erccs, k = 1)
assay(umi.qc, "ruvg1") <- log2(
  t(t(ruvg$normalizedCounts) / colSums(ruvg$normalizedCounts) * 1e6) + 1
)
ruvg <- RUVg(counts(umi.qc), erccs, k = 10)
assay(umi.qc, "ruvg10") <- log2(
  t(t(ruvg$normalizedCounts) / colSums(ruvg$normalizedCounts) * 1e6) + 1
)


scIdx <- matrix(-1, ncol = max(table(umi.qc$individual)), nrow = 3)
tmp <- which(umi.qc$individual == "NA19098")
scIdx[1, 1:length(tmp)] <- tmp
tmp <- which(umi.qc$individual == "NA19101")
scIdx[2, 1:length(tmp)] <- tmp
tmp <- which(umi.qc$individual == "NA19239")
scIdx[3, 1:length(tmp)] <- tmp
cIdx <- rownames(umi.qc)
ruvs <- RUVs(counts(umi.qc), cIdx, k = 1, scIdx = scIdx, isLog = FALSE)
assay(umi.qc, "ruvs1") <- log2(
  t(t(ruvs$normalizedCounts) / colSums(ruvs$normalizedCounts) * 1e6) + 1
)
ruvs <- RUVs(counts(umi.qc), cIdx, k = 10, scIdx = scIdx, isLog = FALSE)
assay(umi.qc, "ruvs10") <- log2(
  t(t(ruvs$normalizedCounts) / colSums(ruvs$normalizedCounts) * 1e6) + 1
)


combat_data <- logcounts(umi.qc)
mod_data <- as.data.frame(t(combat_data))
# Basic batch removal
mod0 = model.matrix(~ 1, data = mod_data) 
# Preserve biological variability
mod1 = model.matrix(~ umi.qc$individual, data = mod_data) 
# adjust for total genes detected
mod2 = model.matrix(~ umi.qc$total_features, data = mod_data)
assay(umi.qc, "combat") <- ComBat(
  dat = t(mod_data), 
  batch = factor(umi.qc$batch), 
  mod = mod0,
  par.prior = TRUE,
  prior.plots = FALSE
)


assay(umi.qc, "combat_tf") <- ComBat(
  dat = t(mod_data), 
  batch = factor(umi.qc$batch), 
  mod = mod2,
  par.prior = TRUE,
  prior.plots = FALSE
)


glm_fun <- function(g, batch, indi) {
  model <- glm(g ~ batch + indi)
  model$coef[1] <- 0 # replace intercept with 0 to preserve reference batch.
  return(model$coef)
}
effects <- apply(
  logcounts(umi.qc), 
  1, 
  glm_fun, 
  batch = umi.qc$batch, 
  indi = umi.qc$individual
)
corrected <- logcounts(umi.qc) - t(effects[as.numeric(factor(umi.qc$batch)), ])
assay(umi.qc, "glm") <- corrected


## BIOLOGICAL ANALYSIS

library(pcaMethods)
library(pcaReduce)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
set.seed(1234567)

# data normalized with cpm, in log counts we have log transform of the data
# cell type 2 annotation of the data
# cell type 1 merge annotation
deng <- readRDS("deng/deng-reads.rds")
deng
table(colData(deng)$cell_type2)

plotPCA(deng, colour_by = "cell_type2")

deng <- sc3_estimate_k(deng)
metadata(deng)$sc3$k_estimation
plotPCA(deng, colour_by = "cell_type1")

deng <- sc3(deng, ks = 10, biology = TRUE)
deng
sc3_plot_consensus(deng, k = 10, show_pdata = "cell_type2")

sc3_plot_silhouette(deng, k = 10)
