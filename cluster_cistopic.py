#!/usr/bin/env python3
import os
import pickle
import argparse
import matplotlib.pyplot as plt
import pandas as pd
from pycisTopic.clust_vis import find_clusters, run_umap, run_tsne, plot_metadata

def main():
    parser = argparse.ArgumentParser(description="Cluster CistopicObject and plot metadata on UMAP/t-SNE")
    parser.add_argument("--cistopic_obj", "-c", required=True, help="Path to CistopicObject pickle file")
    parser.add_argument("--out_dir", "-o", required=True, help="Directory to save plots and updated object")
    args = parser.parse_args()

    # --- Step 1: Load CistopicObject ---
    with open(args.cistopic_obj, "rb") as f:
        cistopic_obj = pickle.load(f)

    os.makedirs(args.out_dir, exist_ok=True)

    # --- Step 2: Ensure cell_data exists and covers all cell_names ---
    if cistopic_obj.cell_data is None:
        cistopic_obj.cell_data = pd.DataFrame(index=cistopic_obj.cell_names)
    else:
        # Add missing cells
        missing_cells = set(cistopic_obj.cell_names) - set(cistopic_obj.cell_data.index)
        if missing_cells:
            missing_df = pd.DataFrame(index=list(missing_cells), columns=cistopic_obj.cell_data.columns)
            cistopic_obj.cell_data = pd.concat([cistopic_obj.cell_data, missing_df])
        # Reorder
        cistopic_obj.cell_data = cistopic_obj.cell_data.loc[cistopic_obj.cell_names]

    # --- Step 3: Extract metadata from object ---
    # sample_id
    if hasattr(cistopic_obj, "cell_data") and "sample_id" in cistopic_obj.cell_data.columns:
        cistopic_obj.cell_data["sample_id_aligned"] = cistopic_obj.cell_data["sample_id"]
    else:
        print("Warning: sample_id not found in cell_data, filling with 'Unknown'")
        cistopic_obj.cell_data["sample_id_aligned"] = "Unknown"

    # celltype_scrna
    if hasattr(cistopic_obj, "cell_data") and "celltype_scrna" in cistopic_obj.cell_data.columns:
        cistopic_obj.cell_data["celltype_scrna_aligned"] = cistopic_obj.cell_data["celltype_scrna"]
    else:
        print("Warning: celltype_scrna not found in cell_data, filling with 'Unknown'")
        cistopic_obj.cell_data["celltype_scrna_aligned"] = "Unknown"

    # --- Step 4: Clustering ---
    print("Running clustering...")
    find_clusters(
        cistopic_obj,
        target='cell',
        k=10,
        res=[0.6, 1.2, 3],
        prefix='pycisTopic_',
        scale=True,
        split_pattern='-'
    )

    # --- Step 5: Compute embeddings ---
    print("Computing UMAP embedding...")
    run_umap(cistopic_obj, target='cell', scale=True)

    print("Computing t-SNE embedding...")
    run_tsne(cistopic_obj, target='cell', scale=True)

    # --- Step 6: Prepare variables to plot ---
    clustering_columns = [col for col in cistopic_obj.cell_data.columns if col.startswith("pycisTopic_")]
    metadata_columns = ["sample_id_aligned", "celltype_scrna_aligned"]
    variables_to_plot = clustering_columns + metadata_columns

    # Filter variables that exist and have at least one non-NaN value
    variables_to_plot = [v for v in variables_to_plot
                         if v in cistopic_obj.cell_data.columns and cistopic_obj.cell_data[v].notna().any()]

    # --- Step 7: Plot metadata on UMAP ---
    if variables_to_plot:
        print(f"Plotting metadata for {len(cistopic_obj.cell_names)} cells: {variables_to_plot}")
        plt.figure(figsize=(12,12))
        plot_metadata(
            cistopic_obj,
            reduction_name='UMAP',
            variables=variables_to_plot,
            target='cell',
            num_columns=4,
            text_size=10,
            dot_size=5
        )
        fig_file = os.path.join(args.out_dir, "metadata_UMAP_all_cells.png")
        plt.savefig(fig_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved metadata UMAP plot: {fig_file}")
    else:
        print("No valid metadata columns found for plotting.")

    # --- Step 8: Save updated CistopicObject ---
    updated_obj_file = os.path.join(args.out_dir, "cistopic_obj_with_clusters_and_metadata.pkl")
    with open(updated_obj_file, "wb") as f:
        pickle.dump(cistopic_obj, f)
    print(f"Saved updated CistopicObject with clustering and metadata: {updated_obj_file}")

    print("All done!")

if __name__ == "__main__":
    main()

