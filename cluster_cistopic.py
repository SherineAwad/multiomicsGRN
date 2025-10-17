#!/usr/bin/env python3
import argparse
import pickle
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np

def cluster_cistopic_manual(input_pickle, output_pickle, outdir, resolutions=[0.6, 1.2, 3], k=10):
    os.makedirs(outdir, exist_ok=True)

    # Load object
    with open(input_pickle, "rb") as f:
        obj = pickle.load(f)

    print(f"Available cell types: {obj.cell_data['celltype'].unique()}")

    # MANUAL CLUSTERING using scanpy
    print("Performing manual clustering with Scanpy...")
    
    # Create AnnData object from cell_topic matrix
    adata = sc.AnnData(X=obj.cell_topic.values)
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    adata.var_names = [f"Topic_{i+1}" for i in range(adata.n_vars)]
    
    # Add metadata
    for col in obj.cell_data.columns:
        adata.obs[col] = obj.cell_data[col].values
    
    # Standard preprocessing
    sc.pp.neighbors(adata, n_neighbors=k, use_rep='X')
    
    # Clustering at different resolutions
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=f'pycisTopic_leiden_{res}')
    
    # UMAP
    sc.tl.umap(adata)
    
    # Transfer results back to cistopic object
    for res in resolutions:
        cluster_col = f'pycisTopic_leiden_{res}'
        obj.cell_data[cluster_col] = adata.obs[cluster_col].values
    
    obj.cell_data['UMAP_1'] = adata.obsm['X_umap'][:, 0]
    obj.cell_data['UMAP_2'] = adata.obsm['X_umap'][:, 1]

    # --- PLOTTING ---
    
    # 1. CELL TYPE UMAP - with cell type labels on plot
    plt.figure(figsize=(12, 10))
    
    unique_celltypes = obj.cell_data['celltype'].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_celltypes)))
    color_map = dict(zip(unique_celltypes, colors))
    
    # Plot each cell type and calculate centroids for labels
    centroids = {}
    for celltype in unique_celltypes:
        mask = obj.cell_data['celltype'] == celltype
        x_data = obj.cell_data.loc[mask, 'UMAP_1']
        y_data = obj.cell_data.loc[mask, 'UMAP_2']
        
        plt.scatter(x_data, y_data, c=[color_map[celltype]], 
                   label=celltype, s=10, alpha=0.7)
        
        # Calculate centroid for label placement
        centroids[celltype] = (x_data.mean(), y_data.mean())
    
    # Add cell type labels at centroids
    for celltype, (x, y) in centroids.items():
        plt.annotate(celltype, (x, y), xytext=(5, 5), 
                    textcoords='offset points', fontsize=9)
    
    plt.title('Cell Types - UMAP')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "celltype_umap.png"), bbox_inches='tight', dpi=300)
    plt.close()
    print("✓ Created celltype UMAP with labels")

    # 2. ANNOTATED CLUSTERS - SHOW CELL TYPES ON CLUSTERS
    for res in resolutions:
        cluster_col = f'pycisTopic_leiden_{res}'
        
        plt.figure(figsize=(12, 10))
        
        unique_clusters = sorted(obj.cell_data[cluster_col].unique())
        cluster_centroids = {}
        
        # Plot clusters colored by cell type composition
        for cluster in unique_clusters:
            mask = obj.cell_data[cluster_col] == cluster
            
            # Get dominant cell type for this cluster
            celltype_counts = obj.cell_data.loc[mask, 'celltype'].value_counts()
            dominant_celltype = celltype_counts.index[0] if len(celltype_counts) > 0 else 'Unknown'
            
            # Use the cell type color
            plt.scatter(obj.cell_data.loc[mask, 'UMAP_1'], 
                       obj.cell_data.loc[mask, 'UMAP_2'], 
                       c=[color_map.get(dominant_celltype, 'gray')], 
                       s=10, alpha=0.7)
            
            # Calculate cluster centroid for label
            x_centroid = obj.cell_data.loc[mask, 'UMAP_1'].mean()
            y_centroid = obj.cell_data.loc[mask, 'UMAP_2'].mean()
            cluster_centroids[cluster] = (x_centroid, y_centroid, dominant_celltype)
        
        # ANNOTATE CLUSTERS WITH CELL TYPE LABELS - NO BOXES
        for cluster, (x, y, celltype) in cluster_centroids.items():
            plt.annotate(f'{celltype}', (x, y), xytext=(5, 5), 
                        textcoords='offset points', fontsize=11, fontweight='bold')
        
        plt.title(f'Clusters annotated with dominant cell types (resolution {res})')
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"annotated_clusters_res_{res}.png"), 
                   bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"✓ Created annotated clusters for resolution {res}")

    # 3. TOPIC-CELLTYPE HEATMAP
    try:
        cell_topic_df = obj.cell_topic.copy()
        cell_topic_df['celltype'] = obj.cell_data['celltype'].values
        topic_by_celltype = cell_topic_df.groupby('celltype').mean()
        
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

    print(f"✓ Manual clustering completed and saved to: {output_pickle}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster CistopicObject with manual method")
    parser.add_argument("-i", "--input_pickle", required=True)
    parser.add_argument("-o", "--output_pickle", required=True)
    parser.add_argument("-d", "--outdir", required=True)
    parser.add_argument("--resolutions", nargs="+", type=float, default=[0.6, 1.2, 3])
    parser.add_argument("--k", type=int, default=10)
    args = parser.parse_args()
    cluster_cistopic_manual(args.input_pickle, args.output_pickle, args.outdir, args.resolutions, args.k)
