import os
import argparse
import pandas as pd
import pycisTopic
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk

print("pycisTopic version:", pycisTopic.__version__)
parser = argparse.ArgumentParser(description="Load cell annotation CSV")
parser.add_argument("input_barcode", help="Path to cell annotation CSV")  # positional argument
args = parser.parse_args()
barcode_mapping = args.input_barcode
# -----------------------------
# 1. Define fragments dictionary
# -----------------------------
fragments_dict = {
    "Control": "TH1_atac_fragments.tsv.gz",
    "KO": "TH2_atac_fragments.tsv.gz"
}

# -----------------------------
# 2. Load chromsizes
# -----------------------------
chromsizes = pd.read_table(
    "https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes",
    header=None,
    names=["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)
print(chromsizes.head())

# -----------------------------
# 3. Create output directories
# -----------------------------
out_dir = os.getcwd()

consensus_dir = os.path.join(out_dir, "consensus_peak_calling")
bed_dir       = os.path.join(consensus_dir, "pseudobulk_bed_files")
bw_dir        = os.path.join(consensus_dir, "pseudobulk_bw_files")

os.makedirs(consensus_dir, exist_ok=True)
os.makedirs(bed_dir, exist_ok=True)
os.makedirs(bw_dir, exist_ok=True)

# -----------------------------
# 4. Load barcode-to-celltype CSV
# -----------------------------
# Expected columns: barcode, celltype
cell_data = pd.read_csv(barcode_mapping)
print("Columns in cell_data:", cell_data.columns)
print(cell_data.head())

# -----------------------------
# 5. Create 'sample' column from barcode prefix
# -----------------------------
# Your barcodes should look like: "Control_AAACGAAAGAGCCTA-1" or "KO_TTGGGAAAGTAGCTC-1"
cell_data['sample'] = cell_data['barcode'].str.split("_").str[0]

# Verify all sample names exist in fragments_dict
missing_samples = set(cell_data['sample'].unique()) - set(fragments_dict.keys())
if missing_samples:
    raise ValueError(f"No fragment file provided for these samples: {missing_samples}")

print("Unique samples:", cell_data['sample'].unique())

# -----------------------------
# 6. Run pseudobulk export
# -----------------------------
bw_paths, bed_paths = export_pseudobulk(
    input_data = cell_data,
    variable = "celltype",       # column with cell type labels
    sample_id_col = "sample",    # column with sample names
    chromsizes = chromsizes,
    bed_path = bed_dir,
    bigwig_path = bw_dir,
    path_to_fragments = fragments_dict,
    n_cpu = 10,
    normalize_bigwig = True,
    temp_dir = "/tmp",
    split_pattern = "-"
)

# -----------------------------
# 7. Save bw_paths and bed_paths for reference
# -----------------------------
with open(os.path.join(consensus_dir, "bw_paths.tsv"), "wt") as f:
    for v in bw_paths:
        f.write(f"{v}\t{bw_paths[v]}\n")

with open(os.path.join(consensus_dir, "bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        f.write(f"{v}\t{bed_paths[v]}\n")

print("Pseudobulk export completed successfully!")

