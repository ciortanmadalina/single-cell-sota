# Single cell RNA seq method comparaison

Run the project by following the next steps:

1. Clone github repo:

git clone https://github.com/ciortanmadalina/single-cell-sota.git

2. Start docker and mount github code:

sudo docker run -d -p 8787:8787 -v $PWD/single-cell-sota:/home/rstudio/single-cell-sota quay.io/hemberg-group/scrna-seq-course-rstudio 

Ths project compares the following methods:
- Seurat
- Monocle
- CIDR https://github.com/VCCRI/CIDR
- SC3 : freezes during clustering phase

on the following datasets:

- deng (https://www.ncbi.nlm.nih.gov/pubmed/24408435 : 
To illustrate clustering of scRNA-seq data, we consider the Deng dataset of cells from developing mouse embryo (Deng et al. 2014). We have preprocessed the dataset and created a SingleCellExperiment object in advance. We have also annotated the cells with the cell types identified in the original publication (it is the cell_type2 column in the colData slot). )

- Human Pancreatic Islet scRNA-Seq Dataset
In this dataset there are 60 cells in 6 cell types after we exclude undefined cells and bulk RNA-Seq samples. Reference for the human pancreatic islet dataset:
Li, J. et al. Single-cell transcriptomes reveal characteristic features of human pancreatic islet cell types. EMBO Reports 17, 178–187 (2016).
- Human Brain scRNA-Seq Dataset
 In this dataset there are 420 cells in 8 cell types after we exclude hybrid cells. Reference for the human brain dataset: Darmanis, S. et al. A survey of human brain transcriptome diversity at the single cell level. Proceedings of the National Academy of Sciences 112, 7285–7290 (2015).

- C-elegans


- Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor.(wget -c http://cf.10xgenomics.com/samples/cell-exp/2.0.1/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz)

