# set directory if needed
setwd("D:/courses/scrnaseq/paper")

# install and load libraries
library(dplyr)
library(Seurat)
library(SingleR)
library(celldex)

# Some of the different data format are given below. Load Different data formats of scrnaseq data. 

# 1. .RDS format
rds_obj <- readRDS('name_of_file.rds')

#--------------------------------------------------------------#

# 2. 10x CellRanger .HD5 format
hdf5_obj <- Read10X_h5(filename = 'name_of_file.h5',
                       use.names = TRUE,
                       unique.features = TRUE)
# Create Seurat object
seurat_hd5_obj <- CreateSeuratObject(counts = hd5_obj)


#--------------------------------------------------------------#

# 3. .mtx file
# matrix: count matrix
# features: genes
# barcodes: cells

mtx_obj <- ReadMTX(mtx = "matrix.mtx.gz",
                   features = "features.tsv.gz",
                   cells = "barcodes.tsv.gz")

seurat_mtx_obj <- CreateSeuratObject(mtx_obj)

#--------------------------------------------------------------#


# Demo scrnaseq dataset: GSE145926	(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4339769)
# Single-cell landscape of bronchoalveolar immune cells in COVID-19 patients
# Following demo is of a single scrnaseq sample from the dataset in hd5 format
# Sample GSM4339769 BALF, C141 (scRNA-seq)

# set directory if needed
setwd("D:/courses/scrnaseq/paper")

# Load libraries
library(dplyr)
library(Seurat)

hdf5_obj <- Read10X_h5(filename = 'GSM4339769_C141_filtered_feature_bc_matrix.h5',
                       use.names = TRUE,
                       unique.features = TRUE)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = hdf5_obj)
seurat_obj

# ---Quality control---

View(seurat_obj@meta.data)

# Calculate mitochondrial percentage in all cells and save in new column
# Mitochondrial genes start with "MT-"
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
View(seurat_obj@meta.data)

# ---Filtering---
# Filtering cells based on genes and mitochondrial percent in a cell
# Filter cells that have genes > 200 & < 2500 & percent.mt < 5
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# ---Normalize data---
seurat_obj <- NormalizeData(seurat_obj)

# ---Identify variable features---
# Identifying 2000 most highly variable genes 
seurat_obj <-FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# ---Scaling the data---
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes)


# ---Perform PCA for linear dimensionality reduction---
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Visualize 
DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)

# Determine dimensionality of the data
ElbowPlot(seurat_obj)

# ---Clustering the data---
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
View(seurat_obj@meta.data)

# Resolution
# The resolution controls cluster granularity.
# Higher resolution: more cluster, finer separation
# Lower resolution: fewer clusters, broder groups
# Visualize the clusters using dimplot function and select the best resolution w.r.t data
# 0.5 is the common default resolution used
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
DimPlot(seurat_obj, group.by = "RNA_snn_res.0.5", label = TRUE)

# ---Perform UMAP for non-linear dimensionality reduction---
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
umap_clusters <-DimPlot(seurat_obj, reduction = "umap")

# View cluster IDs
table(Idents(seurat_obj))

# ---Find differentially expressed features/markers---

# Find identity of a cluster (find which cell type it represents by identifying the markers of a specific cluster.: findAllMarkers()
# Function: FindAllMarkers() -> This function identifies differentially expressed genes in a cluster and compares it with other cluster

all_markers <- FindAllMarkers(seurat_obj,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)
View(all_markers)

# top markers per cluster
top10 <- all_markers %>% 
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 10)

View(top10)

### Eg. Common markers ###
# T cells → CD3D, CD3E
# CD4 T cells → IL7R, LTB
# CD8 T cells → NKG7, CCL5
# B cells → MS4A1, CD79A
# NK cells → NKG7, GNLY
# Monocytes → LST1, S100A8
# Macrophages → C1QA, APOE
# Dendritic cells → FCER1A, CST3
# Plasma cells → MZB1, JCHAIN
# Proliferating cells → MKI67, TOP2A
# Epithelial cells → EPCAM, KRT18

# Visualize marker genes
FeaturePlot(seurat_obj, 
            features = c("CD3D", "MS4A1", "NKG7"))
VlnPlot(seurat_obj, 
        features = c("CD3D", "MS4A1", "NKG7"))

# If the known marker are coming up in the top10 markerr of the clusters, label the cluster with corresponding cell type.

# Another way is matchin with reference databases
# SingleR package

# Extract normalized data
expr <- GetAssayData(seurat_obj, slot = "data")

# Load reference data (Human Primary cell atlas data directly available in the celldex package)
ref <- HumanPrimaryCellAtlasData()

# run automatic cell type annotation
predictions <- SingleR(
  test = expr,
  ref = ref,
  labels = ref$label.main
)

# Add annotations to seurat metadata
seurat_obj$celltype <- predictions$labels


# Visualize annotated clusters
annotated_clusters <- DimPlot(
  seurat_obj,
  group.by = "celltype",
  label = TRUE
)

# Visualize umap and annotated clusters

umap_clusters + annotated_clusters

# ---Save the seurat object---

saveRDS(seurat_obj, "seurat_demo_object.rds")




