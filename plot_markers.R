library(Signac)
library(Seurat)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)
library(argparse)
library(SeuratDisk)

parser <- ArgumentParser(description = "Process Seurat/ATAC data")
parser$add_argument("markers", help = "Text file with one marker gene per line")

args <- parser$parse_args()
markers <- args$markers
gene_markers <- readLines(markers)

gene_markers

myObject <- readRDS("Neurog2_ATAC.rds") 
gene.activities <- GeneActivity(myObject)
myObject[["RNA"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(myObject) <- "RNA"
myObject <- NormalizeData(myObject)
myObject <- ScaleData(myObject)

valid_genes <- gene_markers[gene_markers %in% rownames(myObject)]

for (gene in valid_genes) {
  if (gene %in% rownames(myObject)) {
    # Set the output filename for each gene
    png_filename <- paste0("ATAC_FeaturePlot_", gene, ".png")

    # Open PNG device
    png(filename = png_filename, width = 6, height = 5, units = "in", res = 300)

    # Create the plot
    p <- FeaturePlot(
      myObject,
      features = gene,
      reduction = "umap"
    ) + ggtitle(gene) + theme(plot.title = element_text(hjust = 0.5))

    # Print to the PNG device
    print(p)

    # Close the PNG device
    dev.off()
  }
}

saveRDS(myObject, file = "Neurog2_ATAC_Genes.rds")




