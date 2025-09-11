# -----------------------------
# Full corrected TH1/TH2 ATAC-seq workflow
# -----------------------------
library(Seurat)
library(Signac)
library(GenomicRanges)
library(Matrix)
library(dplyr)
library(plyr)
library(EnsDb.Mmusculus.v79)
library(future)
library(argparse)

# ----------------------------------------
# 1. Load Seurat object
# ----------------------------------------
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assume the first argument is the path to the RDS file
object_name <- args[1]

# -----------------------------
# 1. Parallelization for speed
# -----------------------------
#plan("multisession", workers = 6)  # adjust cores per your SLURM job
plan(sequential)  # prevents multiple copies in RAM

# Increase max globals size for large Seurat objects
options(future.globals.maxSize = 20 * 1024^3)  # 20 GB

# -----------------------------
# 2. Define files
# -----------------------------
frags_control_file <- "TH1_atac_fragments.tsv.gz"
frags_ko_file      <- "TH2_atac_fragments.tsv.gz"

matrix_control_file <- "TH1_filtered_feature_bc_matrix.h5"
matrix_ko_file      <- "TH2_filtered_feature_bc_matrix.h5"

# -----------------------------
# 3. Read 10x ATAC peak matrices
# -----------------------------
counts_control <- Read10X_h5(matrix_control_file)
counts_ko      <- Read10X_h5(matrix_ko_file)

if (is.list(counts_control)) counts_control <- counts_control$Peaks
if (is.list(counts_ko))      counts_ko      <- counts_ko$Peaks

# Fix rownames to Signac format (chr-start-end)
rownames(counts_control) <- gsub(":", "-", rownames(counts_control))
rownames(counts_ko)      <- gsub(":", "-", rownames(counts_ko))

# -----------------------------
# 4. Create ChromatinAssay objects
# -----------------------------
control_assay <- CreateChromatinAssay(
  counts = counts_control,
  fragments = frags_control_file,
  genome = "mm10"
)

ko_assay <- CreateChromatinAssay(
  counts = counts_ko,
  fragments = frags_ko_file,
  genome = "mm10"
)

# -----------------------------
# 5. Create Seurat objects
# -----------------------------
control_seurat <- CreateSeuratObject(counts = control_assay, assay = "peaks")
ko_seurat      <- CreateSeuratObject(counts = ko_assay, assay = "peaks")

control_seurat$condition <- "Control"
ko_seurat$condition      <- "KO"

# -----------------------------
# 6. Merge Seurat objects
# -----------------------------
myObject <- merge(
  control_seurat,
  y = ko_seurat,
  add.cell.ids = c("Control", "KO")
)

# -----------------------------
# 7. Add gene annotations for TSS enrichment
# -----------------------------
annotations <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79)

# Fix chromosome names and genome to match peaks
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
genome(annotations) <- "mm10"

Annotation(myObject) <- annotations

# -----------------------------
# 8. QC metrics
# -----------------------------
myObject <- NucleosomeSignal(myObject)

# TSSEnrichment (will now work with large objects)
myObject <- TSSEnrichment(myObject)

myObject <- subset(
  myObject,
  subset = nCount_peaks < 50000 &
           nucleosome_signal < 4 &
           TSS.enrichment > 2
)

# -----------------------------
# 9. TF-IDF normalization, LSI, clustering
# -----------------------------
myObject <- RunTFIDF(myObject)
myObject <- FindTopFeatures(myObject)
myObject <- RunSVD(myObject)

myObject <- FindNeighbors(myObject, reduction = 'lsi', dims = 2:30)
myObject <- FindClusters(myObject, resolution = 0.5)
myObject <- RunUMAP(myObject, reduction = 'lsi', dims = 2:30)

# ----------------------------------------
# 10. Visualization
# ----------------------------------------
# Define output file
png(filename = "ATAC_clusters.png", width = 2000, height = 1600, res = 300)
# Generate the plot
DimPlot(myObject,
        reduction = "umap",
        group.by = "seurat_clusters",
        label = TRUE,
        repel = TRUE) +
  ggtitle("ATAC")

# Close the device to save the file
dev.off()
png(filename = "ATAC_samples.png", width = 2000, height = 1600, res = 300)
DimPlot(
  myObject,
  reduction = "umap",
  group.by = "condition",
  cols = c("Control" = "blue", "KO" = "red"),
  label = TRUE,
  repel = TRUE
)
dev.off()

# -----------------------------
# 11. Save Seurat object
# -----------------------------
saveRDS(myObject, file = object_name) 

