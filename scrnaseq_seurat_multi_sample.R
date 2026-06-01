library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

#Read all files,add sample ids and create seurat objects

seurat_h5_list <- function(dir_path){
  files <- list.files(path = dir_path, pattern = "*.h5", full.names = TRUE)
  
  sample_names <- tools::file_path_sans_ext(basename(files))
  
  seurat_list <- lapply(seq_along(files), function(i){
  
     data <- Read10X_h5(files[i])
   
     obj <- CreateSeuratObject(
       counts = data,
       project = sample_names[i]
   )
     
     obj$sample <- sample_names[i]

      return(obj)
  })
  
  names(seurat_list) <- sample_names
  
  return(seurat_list)
  
}
############################################################################
# Call function
seurat_list <- seurat_h5_list("D:/courses/scrnaseq/paper/multisample")
seurat_list
before_qc_cells = sapply(seurat_list, ncol)
for (sample in names(before_qc_cells)){
  cat(
    "Sample: ", sample,
    "| Cells before QC: ", before_qc_cells[sample],
    "\n"
  )
}
#write.csv(before_qc_cells, "before_qc_cells.csv")
############################################################################

# Calculate mitochondrial percent and  QC

qc_seurat <- function(obj, min_features, max_features, max_mt){
  
  #calculate mitochondrial percentage
  obj[["percent.mt"]] <- PercentageFeatureSet(obj,pattern = "^MT-" )
  
  # Filtering
  obj <- subset(
    obj, subset = nFeature_RNA >= min_features &
                  nFeature_RNA <= max_features &
                  percent.mt <= max_mt
  )
  
  return(obj)
} 

############################################################################
#call function
seurat_list_qc <- lapply(seurat_list, qc_seurat, min_features = 200, max_features = 2500, max_mt = 5)
seurat_list_qc
after_qc_cells <- sapply(seurat_list_qc, ncol)

for (sample in names(after_qc_cells)){
  cat(
    "Sample: ", sample,
    "| Cells after QC: ", after_qc_cells[sample],
    "\n"
  )
}

############################################################################
#Merging objects and preprocessing

names(seurat_list_qc)

merged_obj <- function(obj, project_name){
  cat("Object:\n")
  print(obj)
  
  cat("\nProject name: ", project_name, "\n")
  
  combined <- merge(
    x = obj[[1]],
    y = obj[-1],
    add.cell.ids = names(obj),
    project = project_name
  )
  
  return(combined)
}

############################################################################
#call function
combined <- merged_obj(obj = seurat_list_qc, project_name = "COVID_BALF")
ncol(combined)
sum(sapply(seurat_list_qc, ncol))
############################################################################

#Visualize combined object clustering before integration

comb1 <- NormalizeData(combined)
comb <- FindVariableFeatures(comb1, selection.method = "vst", nfeatures = 2000)
comb <- ScaleData(comb)
comb <- RunPCA(comb, npcs = 30)
ElbowPlot(comb, ndims = 30, reduction = "pca")
comb <- RunUMAP(comb, dims = 1:30)
plot1 <- DimPlot(comb, group.by="orig.ident")
plot1

############################################################################

# Integration using seurat

# Normalization
seurat_list_qc <- lapply(seurat_list_qc, NormalizeData)

# Find Variable features
seurat_list_qc <- lapply(seurat_list_qc, 
                         FindVariableFeatures, 
                         selection.method = "vst",
                         nfeatures = 2000
)

# Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_list_qc, nfeatures = 3000)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list_qc, anchor.features = features)

#Integration
integrated <- IntegrateData(anchorset = anchors)

#Scaling
integrated <- ScaleData(integrated)

#PCA
integrated <- RunPCA(integrated)

#UMAP
integrated <- RunUMAP(integrated, dims = 1:30)

plot2 <- DimPlot(integrated, group.by="orig.ident")
plot2
plot1+plot2
