# Load first dataset sce10x_qc

load("single-cell-sota/input/cellBench/sincell_with_class.RData")

counts(sce10x_qc)[1:5, 1:5]
head(colData(sce10x_qc))
table(colData(sce10x_qc)$cell_line)

write.csv(counts(sce10x_qc),'single-cell-sota/input/cellBench/sce10x_qc.csv')
write.csv(colData(sce10x_qc)$cell_line,'single-cell-sota/input/cellBench/sce10x_qc_truth.csv')


rm(list = ls())
.rs.restartR()


load("single-cell-sota/input/cellBench/mRNAmix_qc.RData")

write.csv(counts(sce8_qc),'single-cell-sota/input/cellBench/sce8_qc.csv')
write.csv(sce8_qc$mix,'single-cell-sota/input/cellBench/sce8_qc_truth.csv')

write.csv(counts(sce2_qc),'single-cell-sota/input/cellBench/sce2_qc.csv')
write.csv(sce2_qc$mix,'single-cell-sota/input/cellBench/sce2_qc_truth.csv')



rm(list = ls())
.rs.restartR()


load("single-cell-sota/input/cellBench/9cellmix_qc.RData")
