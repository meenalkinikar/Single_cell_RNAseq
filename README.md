# Single_Cell_RNAseq

A reproducible single-cell RNA-seq analysis workflow developed in R using the Seurat package for quality control, normalization, dimensionality reduction, clustering, cell-type identification, and downstream transcriptomic analysis. The pipeline includes filtering low-quality cells, removal of potential doublets and high-mitochondrial-content cells, data normalization, identification of highly variable genes, principal component analysis (PCA), graph-based clustering, UMAP/t-SNE visualization, marker gene identification, and differential gene expression analysis to characterize cellular heterogeneity within complex biological samples.

## Workflow

1. Data import and quality control
2. Filtering low-quality cells and genes
3. Mitochondrial content assessment
4. Data normalization and scaling
5. Identification of highly variable genes
6. Principal Component Analysis (PCA)
7. Cell clustering using graph-based methods
8. UMAP/t-SNE dimensionality reduction and visualization
9. Cluster marker gene identification
10. Cell type annotation using canonical markers
11. Differential gene expression analysis
12. Visualization of cellular heterogeneity and population structure

## Please Note

Sample data has been provided for demonstration purposes. The analysis scripts can be executed using R and the Seurat package. The workflow is designed to generate publication-quality visualizations and identify biologically meaningful cell populations from single-cell transcriptomic datasets.

## Multisample analysis
Data download link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926&utm_source=chatgpt.com





