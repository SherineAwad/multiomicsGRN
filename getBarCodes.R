library(Seurat)
library(argparse)

# 1. Load your clustered ATAC Seurat object
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assume the first argument is the path to the RDS file
object_path <- args[1]
atac_barcode_path <- args[2]

# Read the RDS file
myObject <- readRDS(object_path)

# 2. Extract barcodes and cluster IDs
# 'Idents' contains the cluster assignments
cell_annotation <- data.frame(
  barcode = colnames(myObject),
  celltype = Idents(myObject)
)

# 3. (Optional) add sample information if you have it
# If your barcodes contain sample prefixes like "Control_" or "KO_", you can extract it:
cell_annotation$sample <- sapply(strsplit(as.character(cell_annotation$barcode), "_"), `[`, 1)

# 4. Save as CSV
write.csv(cell_annotation,atac_barcode_path, row.names = FALSE)

# 5. Verify
head(cell_annotation)

