#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import argparse
import os

def main(input_file, expr_out, gene_meta_out, cell_meta_out):
    # Load h5ad
    print(f"Loading {input_file}...")
    adata = sc.read_h5ad(input_file)

    # -------------------------------
    # 1. Gene expression matrix
    # Use existing normalized/log-transformed data
    print("Extracting gene expression matrix...")
    if hasattr(adata.X, "toarray"):
        expr_matrix = pd.DataFrame(adata.X.toarray(),
                                   index=adata.obs_names,
                                   columns=adata.var_names)
    else:
        expr_matrix = pd.DataFrame(adata.X,
                                   index=adata.obs_names,
                                   columns=adata.var_names)

    # Save gene expression matrix
    expr_dir = os.path.dirname(expr_out)
    if expr_dir and not os.path.exists(expr_dir):
        os.makedirs(expr_dir)
    expr_matrix.to_csv(expr_out)
    print(f"Saved gene expression matrix to {expr_out}")

    # -------------------------------
    # 2. Gene metadata
    print("Extracting gene metadata...")
    gene_meta = adata.var.copy()  # includes var_names and other gene info
    gene_meta_dir = os.path.dirname(gene_meta_out)
    if gene_meta_dir and not os.path.exists(gene_meta_dir):
        os.makedirs(gene_meta_dir)
    gene_meta.to_csv(gene_meta_out)
    print(f"Saved gene metadata to {gene_meta_out}")

    # -------------------------------
    # 3. Cell metadata
    print("Extracting cell metadata...")
    cell_meta = adata.obs.copy()  # includes clustering and annotations
    cell_meta_dir = os.path.dirname(cell_meta_out)
    if cell_meta_dir and not os.path.exists(cell_meta_dir):
        os.makedirs(cell_meta_dir)
    cell_meta.to_csv(cell_meta_out)
    print(f"Saved cell metadata to {cell_meta_out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate SCENIC+ input files from h5ad")
    parser.add_argument("-i", "--input", required=True, help="Input h5ad file")
    parser.add_argument("--expr_out", required=True, help="Output CSV for gene expression matrix (cells x genes)")
    parser.add_argument("--gene_meta_out", required=True, help="Output CSV for gene metadata")
    parser.add_argument("--cell_meta_out", required=True, help="Output CSV for cell metadata")
    args = parser.parse_args()

    main(args.input, args.expr_out, args.gene_meta_out, args.cell_meta_out)

