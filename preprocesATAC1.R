# ---------------------------
# Load libraries
# ---------------------------
library(Seurat)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomeInfoDb)
library(SeuratDisk)

# ---------------------------
# Step 0: Fragment files
# ---------------------------
frags_control_file <- "TH1_atac_fragments.tsv.gz"  # Control
frags_ko_file      <- "TH2_atac_fragments.tsv.gz"  # KO

# Wrap fragment files as Signac Fragment objects
frags_control <- CreateFragmentObject(frags_control_file)
frags_ko      <- CreateFragmentObject(frags_ko_file)

# ---------------------------
# Step 1: Create genome bins (5kb)
# ---------------------------
genome <- BSgenome.Mmusculus.UCSC.mm10
bins <- tileGenome(
  seqlengths = seqlengths(genome),
  tilewidth = 5000,            # 5kb bins
  cut.last.tile.in.chrom = TRUE
)

# ---------------------------
# Step 2: Create FeatureMatrix (cells x bins)
# ---------------------------
counts_control <- FeatureMatrix(
  fragments = frags_control,
  features  = bins
)

counts_ko <- FeatureMatrix(
  fragments = frags_ko,
  features  = bins
)

# ---------------------------
# Step 3: Create ChromatinAssays
# ---------------------------
assay_control <- CreateChromatinAssay(
  counts    = counts_control,
  fragments = frags_control,
  assay     = "peaks"
)

assay_ko <- CreateChromatinAssay(
  counts    = counts_ko,
  fragments = frags_ko,
  assay     = "peaks"
)

# ---------------------------
# Step 4: Create Seurat objects
# ---------------------------
# Use counts from ChromatinAssay to create Seurat object
seurat_control <- CreateSeuratObject(
  counts  = assay_control@counts,
  assay   = "peaks",
  project = "Neurog2_ATAC"
)
# Attach the full ChromatinAssay
seurat_control[["peaks"]] <- assay_control

seurat_ko <- CreateSeuratObject(
  counts  = assay_ko@counts,
  assay   = "peaks",
  project = "Neurog2_ATAC"
)
seurat_ko[["peaks"]] <- assay_ko

# Optional: Filter low-quality cells
seurat_control <- subset(seurat_control, subset = nCount_peaks > 1000)
seurat_ko      <- subset(seurat_ko, subset = nCount_peaks > 1000)

# ---------------------------
# Step 5: Merge samples
# ---------------------------
combined <- merge(
  x            = seurat_control,
  y            = seurat_ko,
  add.cell.ids = c("Control","KO"),
  project      = "Neurog2_ATAC"
)

# ---------------------------
# Step 6: Filter low-quality bins
# ---------------------------
bin_counts <- rowSums(combined@assays$peaks@counts > 0)
combined <- subset(combined, features = names(bin_counts[bin_counts >= 10]))

# ---------------------------
# Step 7: Binarize counts for SCENIC+
# ---------------------------
combined@assays$peaks@counts <- (combined@assays$peaks@counts > 0) * 1

# ---------------------------
# Step 8: Save
# ---------------------------
SaveH5Seurat(combined, filename = "Neurog2_ATAC.h5Seurat")

