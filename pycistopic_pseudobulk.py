import os
import pandas as pd
import pycisTopic
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk

print("pycisTopic version:", pycisTopic.__version__)

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
print("Chromsizes head:")
print(chromsizes.head())

# -----------------------------
# 3. Create output directories
# -----------------------------
out_dir = "scenicOuts"
consensus_dir = os.path.join(out_dir, "consensus_peak_calling")
bed_dir       = os.path.join(consensus_dir, "pseudobulk_bed_files")
bw_dir        = os.path.join(consensus_dir, "pseudobulk_bw_files")

os.makedirs(consensus_dir, exist_ok=True)
os.makedirs(bed_dir, exist_ok=True)
os.makedirs(bw_dir, exist_ok=True)

# -----------------------------
# 4. Create cell data from fragments files (SKIPPING HEADER LINES)
# -----------------------------
cell_data_list = []

for sample, fragments_file in fragments_dict.items():
    print(f"Reading barcodes from {fragments_file}...")
    
    if not os.path.exists(fragments_file):
        print(f"  ERROR: File not found!")
        continue
    
    try:
        # Read fragments file - SKIP HEADER LINES starting with #
        fragments_df = pd.read_csv(
            fragments_file, 
            sep='\t', 
            header=None,
            comment='#',  # THIS IS THE KEY FIX: skip lines starting with #
            names=['chrom', 'start', 'end', 'barcode', 'count'],
            usecols=[0, 1, 2, 3, 4]
        )
        
        print(f"  Data shape (after skipping headers): {fragments_df.shape}")
        print(f"  First few barcodes: {fragments_df['barcode'].head().tolist()}")
        
        # Get unique barcodes
        barcodes = fragments_df['barcode'].unique()
        print(f"  Found {len(barcodes)} unique barcodes")
        
        # Create DataFrame for this sample
        sample_df = pd.DataFrame({
            'barcode': barcodes.astype(str),
            'sample': str(sample),
            'celltype': str(sample)
        })
        
        cell_data_list.append(sample_df)
        print(f"  Successfully processed {len(barcodes)} barcodes for sample {sample}")
        
    except Exception as e:
        print(f"  Error processing {fragments_file}: {e}")
        # Try alternative approach if the first fails
        try:
            print("  Trying alternative reading method...")
            fragments_df = pd.read_csv(
                fragments_file, 
                sep='\t', 
                header=None,
                comment='#',
                usecols=[3]  # Only read barcode column
            )
            fragments_df.columns = ['barcode']
            barcodes = fragments_df['barcode'].unique()
            print(f"  Alternative method found {len(barcodes)} barcodes")
            
            sample_df = pd.DataFrame({
                'barcode': barcodes.astype(str),
                'sample': str(sample),
                'celltype': str(sample)
            })
            cell_data_list.append(sample_df)
            
        except Exception as e2:
            print(f"  Alternative method also failed: {e2}")

# Check if we have any data
if len(cell_data_list) == 0:
    print("\nERROR: No barcode data found in any fragment file!")
    exit(1)

# Combine all samples
cell_data = pd.concat(cell_data_list, ignore_index=True)

# Create simplified version with only required columns
cell_data_simple = cell_data[['barcode', 'sample', 'celltype']].copy()

# Ensure all columns are strings
for col in ['barcode', 'sample', 'celltype']:
    cell_data_simple[col] = cell_data_simple[col].astype(str)

print(f"\nFinal cell data summary:")
print(f"Total barcodes across all samples: {len(cell_data_simple)}")
print("Sample distribution:")
print(cell_data_simple['sample'].value_counts())
print("\nFirst few rows:")
print(cell_data_simple.head())

# -----------------------------
# 5. Run pseudobulk export
# -----------------------------
print("\nStarting pseudobulk export...")

try:
    bw_paths, bed_paths = export_pseudobulk(
        input_data=cell_data_simple,
        variable="celltype",
        sample_id_col="sample", 
        chromsizes=chromsizes,
        bed_path=bed_dir,
        bigwig_path=bw_dir,
        path_to_fragments=fragments_dict,
        n_cpu=4,
        normalize_bigwig=True,
        temp_dir="/tmp"
    )
    
    print("Pseudobulk export completed successfully!")
    
    # Save paths
    with open(os.path.join(consensus_dir, "bw_paths.tsv"), "w") as f:
        for v in bw_paths:
            f.write(f"{v}\t{bw_paths[v]}\n")

    with open(os.path.join(consensus_dir, "bed_paths.tsv"), "w") as f:
        for v in bed_paths:
            f.write(f"{v}\t{bed_paths[v]}\n")

    print(f"\nOutput directories:")
    print(f"  - BED files: {bed_dir}")
    print(f"  - BigWig files: {bw_dir}")
    print(f"  - Path lists: {consensus_dir}")

    # Print summary of generated files
    print("\nGenerated files summary:")
    for celltype_sample, bed_path in bed_paths.items():
        if os.path.exists(bed_path):
            file_size = os.path.getsize(bed_path) / (1024 * 1024)
            print(f"  - {celltype_sample}: {bed_path} ({file_size:.2f} MB)")
            
except Exception as e:
    print(f"Error during pseudobulk export: {e}")
    
    # Final attempt with minimal parameters
    print("\nTrying with minimal parameters...")
    try:
        bw_paths, bed_paths = export_pseudobulk(
            input_data=cell_data_simple,
            variable="celltype",
            sample_id_col="sample",
            chromsizes=chromsizes,
            bed_path=bed_dir,
            bigwig_path=bw_dir,
            path_to_fragments=fragments_dict,
            n_cpu=2
        )
        print("Success with minimal parameters!")
        
        # Save paths
        with open(os.path.join(consensus_dir, "bw_paths.tsv"), "w") as f:
            for v in bw_paths:
                f.write(f"{v}\t{bw_paths[v]}\n")

        with open(os.path.join(consensus_dir, "bed_paths.tsv"), "w") as f:
            for v in bed_paths:
                f.write(f"{v}\t{bed_paths[v]}\n")
                
    except Exception as e2:
        print(f"Final error: {e2}")
