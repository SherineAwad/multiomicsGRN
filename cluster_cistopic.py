#!/usr/bin/env python3
import argparse
import pickle
import os
import pandas as pd
import seaborn as sns
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    cell_topic_heatmap
)
import matplotlib.pyplot as plt

def cluster_cistopic(input_pickle, output_pickle, outdir, resolutions=[0.6, 1.2, 3], k=10):
    os.makedirs(outdir, exist_ok=True)

    # Load object
    with open(input_pickle, "rb") as f:
        obj = pickle.load(f)

    # PRESERVE original cell_data
    original_cell_data = obj.cell_data.copy()

    # --- Clustering ---
    find_clusters(obj, target='cell', k=k, res=resolutions, prefix='pycisTopic_', scale=True)

    # --- Dimensionality Reduction ---
    run_umap(obj, target='cell', scale=True)
    run_tsne(obj, target='cell', scale=True)

    # FIX: Extract clustering results before restoring
    clustering_cols = [col for col in obj.cell_data.columns if col.startswith('pycisTopic_')]
    clustering_data = obj.cell_data[clustering_cols].copy()

    # Restore original cell_data with correct barcode index
    obj.cell_data = original_cell_data
    obj.cell_data.index = obj.cell_names  # Set barcodes as index

    # Add clustering results back
    for col in clustering_cols:
        obj.cell_data[col] = clustering_data[col].values

    # --- Simple plotting that works ---
    # Plot cell types
    try:
        plot_metadata(obj, reduction_name='UMAP', variables=['celltype_scrna'], target='cell')
        plt.savefig(os.path.join(outdir, "celltype_umap.png"), bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"Cell type plot failed: {e}")

    # REPLACED: Topic UMAP with Topic-Celltype Heatmap
    try:
        # Get cell_topic matrix and cell types
        cell_topic_df = obj.selected_model.cell_topic.T  # Cells x Topics
        cell_topic_df['celltype'] = obj.cell_data['celltype_scrna']
        
        # Group by cell type and average topic contributions
        topic_by_celltype = cell_topic_df.groupby('celltype').mean()
        
        # Plot heatmap
        plt.figure(figsize=(12, 8))
        sns.heatmap(topic_by_celltype, cmap='viridis', annot=False, linewidths=0.5)
        plt.title('Topic Enrichment by Cell Type')
        plt.xlabel('Topics')
        plt.ylabel('Cell Types')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "topic_celltype_heatmap.png"), dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Created topic-celltype heatmap")
    except Exception as e:
        print(f"Topic heatmap failed: {e}")

    # Save object
    with open(output_pickle, "wb") as f:
        pickle.dump(obj, f)

    print(f"✓ Clustering completed and saved to: {output_pickle}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster CistopicObject")
    parser.add_argument("-i", "--input_pickle", required=True)
    parser.add_argument("-o", "--output_pickle", required=True)
    parser.add_argument("-d", "--outdir", required=True)
    parser.add_argument("--resolutions", nargs="+", type=float, default=[0.6, 1.2, 3])
    parser.add_argument("--k", type=int, default=10)
    args = parser.parse_args()
    cluster_cistopic(args.input_pickle, args.output_pickle, args.outdir, args.resolutions, args.k)
