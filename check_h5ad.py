import scanpy as sc
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('h5ad_file', help='Input h5ad file to inspect')
args = parser.parse_args()

# Load the file
adata = sc.read_h5ad(args.h5ad_file)

print("="*50)
print("BASIC INFO:")
print(f"Shape: {adata.shape}")
print(f"Observations (cells): {adata.n_obs}")
print(f"Variables (genes): {adata.n_vars}")
print()

print("OBS (cell metadata):")
print(adata.obs)
print()

# 🧠 Check if the file is annotated with cell type info
celltype_cols = [col for col in adata.obs.columns if "celltype" in col.lower() or "cell_type" in col.lower()]
if len(celltype_cols) > 0:
    print("Detected potential cell type annotation columns:", celltype_cols)
    # Try to detect sample column
    sample_cols = [col for col in adata.obs.columns if "sample" in col.lower()]
    sample_col = sample_cols[0] if len(sample_cols) > 0 else None
    if sample_col:
        print(f"Detected sample column: '{sample_col}'")
    
    for col in celltype_cols:
        print(f"\nCell counts per type in '{col}':")
        print(adata.obs[col].value_counts())
        
        if sample_col:
            print(f"\nCell counts per type in '{col}' per sample:")
            counts_per_sample = adata.obs.groupby([sample_col, col]).size()
            print(counts_per_sample)
else:
    print("No cell type annotation column detected in adata.obs.")
print()

print("VAR (gene metadata):")
print(adata.var)
print()

print("OBSM keys:", list(adata.obsm.keys()))
for key in adata.obsm.keys():
    print(f"  {key}: {adata.obsm[key].shape}")

print("VARSM keys:", list(adata.varm.keys()))
for key in adata.varm.keys():
    print(f"  {key}: {adata.varm[key].shape}")

print("LAYERS keys:", list(adata.layers.keys()))
for key in adata.layers.keys():
    print(f"  {key}: {adata.layers[key].shape}")

print("UNS keys:", list(adata.uns.keys()))
for key in adata.uns.keys():
    print(f"  {key}: {type(adata.uns[key])}")

print("RAW attribute exists:", adata.raw is not None)
if adata.raw is not None:
    print("RAW shape:", adata.raw.shape)

print("="*50)

