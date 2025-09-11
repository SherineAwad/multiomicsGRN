# ----------------------------------------
# Load libraries
# ----------------------------------------
library(Seurat)
library(Signac)
library(SeuratDisk)
library(ggplot2)
library(argparse)

# ----------------------------------------
# 1. Load Seurat object
# ----------------------------------------
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assume the first argument is the path to the RDS file
object_path <- args[1]

# Read the RDS file
myObject <- readRDS(object_path)

# ----------------------------------------
# 2. Dimensionality reduction & clustering
# ----------------------------------------
# TF-IDF normalization
myObject <- RunTFIDF(myObject)

# Select top features
myObject <- FindTopFeatures(myObject, min.cutoff = 'q0')

# Latent Semantic Indexing (LSI)
myObject <- RunSVD(myObject)

# Build neighborhood graph and clusters
myObject <- FindNeighbors(myObject, reduction = "lsi", dims = 2:30)
myObject <- FindClusters(myObject, resolution = 0.5)

# UMAP for visualization
myObject <- RunUMAP(myObject, reduction = "lsi", dims = 2:30)

# ----------------------------------------
# 3. Visualization
# ----------------------------------------
p1 <- DimPlot(myObject, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Clusters by Seurat")
print(p1)

p2 <- DimPlot(myObject, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Control vs KO")
print(p2)

# ----------------------------------------
# 4. (Optional) Save annotated Seurat object
# ----------------------------------------
newObject <- paste0("Clustered", myObject)
saveRDS(myObject, file = newObject)

