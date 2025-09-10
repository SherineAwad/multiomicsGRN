# ----------------------------------------
# Load libraries
# ----------------------------------------
library(Seurat)
library(Signac)
library(SeuratDisk)
library(ggplot2)

# ----------------------------------------
# 1. Load Seurat object
# ----------------------------------------
combined <- LoadH5Seurat("Neurog2_ATAC.h5Seurat")

# ----------------------------------------
# 2. Dimensionality reduction & clustering
# ----------------------------------------
# TF-IDF normalization
combined <- RunTFIDF(combined)

# Select top features
combined <- FindTopFeatures(combined, min.cutoff = 'q0')

# Latent Semantic Indexing (LSI)
combined <- RunSVD(combined)

# Build neighborhood graph and clusters
combined <- FindNeighbors(combined, reduction = "lsi", dims = 2:30)
combined <- FindClusters(combined, resolution = 0.5)

# UMAP for visualization
combined <- RunUMAP(combined, reduction = "lsi", dims = 2:30)

# ----------------------------------------
# 3. Visualization
# ----------------------------------------
p1 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Clusters by Seurat")
print(p1)

p2 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Control vs KO")
print(p2)

# ----------------------------------------
# 4. Export barcode-to-cluster annotation
# ----------------------------------------
cell_annotation <- data.frame(
  barcode   = colnames(combined),
  cluster   = Idents(combined)  # cluster assignment
)

# Save as CSV (for pycisTopic / SCENIC+ input)
write.csv(cell_annotation, "barcode_to_cluster.csv", row.names = FALSE)

# ----------------------------------------
# 5. (Optional) Save annotated Seurat object
# ----------------------------------------
saveRDS(combined, file = "Neurog2_ATAC_clustered.rds")

