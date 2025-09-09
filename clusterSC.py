import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import os

# Fix for importlib metadata
sys.modules['importlib.metadata'] = importlib_metadata

# -----------------------------
# Parse arguments
# -----------------------------
parser = argparse.ArgumentParser()
parser.add_argument('myObject', help="Input AnnData file (.h5ad)")
parser.add_argument('markers', help="Text file with one marker gene per line")
args = parser.parse_args()

myObject = args.myObject
newObject = "clustered_" + myObject
markers = args.markers
print("Markers file:", markers)

# Base name for figure outputs
base_name = os.path.splitext(os.path.basename(newObject))[0]

# -----------------------------
# Load AnnData
# -----------------------------
combined_adata = sc.read(myObject)

# -----------------------------
# Preprocessing & clustering
# -----------------------------
# If already normalized/logged, you can skip these two lines
sc.pp.normalize_total(combined_adata, target_sum=1e4)
sc.pp.log1p(combined_adata)

# PCA → neighbors → clustering → UMAP
sc.tl.pca(combined_adata, svd_solver="arpack")
sc.pp.neighbors(combined_adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(combined_adata)

# Leiden clustering (with recommended args)
sc.tl.leiden(combined_adata, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)

# Save cluster plot
figure_name = f"_{base_name}_Clusters.png"
sc.pl.umap(combined_adata, color=["leiden"], legend_loc="on data", save=figure_name)

# -----------------------------
# Marker gene plotting
# -----------------------------
with open(markers) as f:
    marker_genes = [line.strip() for line in f]

for gene in marker_genes:
    if gene in combined_adata.var_names:
        sc.pl.umap(
            combined_adata,
            color=gene,
            title=gene,
            save=f"_{base_name}_{gene}.png"
        )
    else:
        print(f"⚠️ Skipping {gene}: not in adata.var_names")

# -----------------------------
# Save clustered object
# -----------------------------
combined_adata.obs_names_make_unique()
combined_adata.write(newObject, compression="gzip")
print(f"✅ Saved clustered AnnData to {newObject}")

