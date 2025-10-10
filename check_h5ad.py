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
