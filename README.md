# Single Cell RNAseq

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

Sample data for single sample analysis has been provided for demonstration purposes. The analysis scripts can be executed using R and the Seurat package. 

################################################################################

## Cell Annotation using Single R

SingleR is an automated cell type annotation method for single-cell RNA sequencing (scRNA-seq) data. It identifies cell types by comparing the gene expression profile of each query cell to annotated reference datasets (e.g., HPCA, Blueprint, Monaco Immune) and assigns the most similar cell type label based on expression similarity.
Pruning in SingleR: After assigning cell types, SingleR removes low-confidence or ambiguous predictions by replacing them with NA in the pruned.labels output. This helps improve annotation reliability by retaining only high-confidence cell type assignments for downstream analysis.

List of reference datasets from celldex:
1. HumanPrimaryCellAtlasData
2. BlueprintEncodeData
3. MonacoImmuneData
4. DatabaseImmuneCellExpressionData (DICE)
5. NovershternHematopoieticData
6. ImmGenData
7. MouseRNAseqData

Best Practice: Treat SingleR annotations as an initial prediction rather than the final annotation. Always validate predicted cell types using canonical marker gene expression, cluster-specific marker analysis (FindAllMarkers()), and, where possible, cross-reference with another annotation method (e.g., Azimuth, CellTypist, or manual annotation) before performing downstream analyses.







