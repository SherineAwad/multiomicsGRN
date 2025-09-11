#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import pycisTopic
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk

print("pycisTopic version:", pycisTopic.__version__)

# -----------------------------
# 1. Parse command-line arguments
# -----------------------------
parser = argparse.ArgumentParser(description="Generate pseudobulk BigWig/BED from ATAC barcode CSV")
parser.add_argument("input_barcode", help="Path to cell annotation CSV")
args = parser.parse_args()
barcode_csv = args.input_barcode

# -----------------------------
# 2. Fragment file paths
# -----------------------------
fragments_dict = {
    "Control": "TH1_atac_fragments.tsv.gz",
    "KO": "TH2_atac_fragments.tsv.gz"
}

# Verify fragment files exist
for sample, path in fragments_dict.items():
    if not os.path.exists(path):
        raise FileNotFoundError(f"Fragment file for {sample} not found: {path}")

# -----------------------------
# 3. Load chromosome sizes (mm10)
# -----------------------------
chromsizes = pd.read_table(
    "https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes",
    header=None,
    names=["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)

# -----------------------------
# 4. Set output directories
# -----------------------------
out_dir = "."
consensus_dir = os.path.join(out_dir, "consensus_peak_calling")
bed_dir       = os.path.join(consensus_dir, "pseudobulk_bed_files")
bw_dir        = os.path.join(consensus_dir, "pseudobulk_bw_files")
temp_dir      = os.path.join(out_dir, "temp_pseudobulk")

# Create directories
os.makedirs(consensus_dir, exist_ok=True)
os.makedirs(bed_dir, exist_ok=True)
os.makedirs(bw_dir, exist_ok=True)
os.makedirs(temp_dir, exist_ok=True)

# -----------------------------
# 5. Load barcode-to-celltype CSV
# -----------------------------
cell_data = pd.read_csv(barcode_csv)

# Fix barcodes: remove "Control_" / "KO_" prefixes so they match fragment files
cell_data['barcode'] = cell_data['barcode'].str.replace(r'^(Control_|KO_)', '', regex=True)

# Convert celltype to string (required by pycisTopic)
cell_data['celltype'] = cell_data['celltype'].astype(str)

# Verify sample names match fragment dictionary keys
missing_samples = set(cell_data['sample'].unique()) - set(fragments_dict.keys())
if missing_samples:
    raise ValueError(f"No fragment file provided for these samples: {missing_samples}")

# -----------------------------
# 6. Run pseudobulk export
# -----------------------------
bw_paths, bed_paths = export_pseudobulk(
    input_data = cell_data,
    variable = "celltype",          # column in CSV for cell type labels
    sample_id_col = "sample",       # column in CSV for sample
    chromsizes = chromsizes,
    bed_path = bed_dir,
    bigwig_path = bw_dir,
    path_to_fragments = fragments_dict,
    n_cpu = 2,
    normalize_bigwig = True,
    temp_dir = temp_dir,
    split_pattern = "-"             # matches barcode suffix "-1"
)

# -----------------------------
# 7. Save BigWig and BED paths
# -----------------------------
with open(os.path.join(consensus_dir, "bw_paths.tsv"), "wt") as f:
    for k in bw_paths:
        f.write(f"{k}\t{bw_paths[k]}\n")

with open(os.path.join(consensus_dir, "bed_paths.tsv"), "wt") as f:
    for k in bed_paths:
        f.write(f"{k}\t{bed_paths[k]}\n")

print("✅ Pseudobulk export completed successfully!")

