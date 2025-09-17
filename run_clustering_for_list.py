#!/usr/bin/env python3
import os
import pickle
import matplotlib.pyplot as plt
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata
)

def main(cistopic_pickle, output_dir):
    # Load MALLET CistopicObject(s)
    with open(cistopic_pickle, "rb") as f:
        cistopic_objs = pickle.load(f)

    # Ensure we have a list
    if not isinstance(cistopic_objs, list):
        cistopic_objs = [cistopic_objs]

    os.makedirs(output_dir, exist_ok=True)

    for idx, obj in enumerate(cistopic_objs):
        sample_name = f"sample_{idx}"
        print(f"Processing {sample_name}...")

        # Use first available categorical metadata if Seurat_cell_type is missing
        metadata_keys = obj.cell_data.columns.tolist()
        if 'Seurat_cell_type' in metadata_keys:
            cell_type_var = 'Seurat_cell_type'
        else:
            # fallback: first categorical column
            cell_type_var = 'sample_id'  

        # 1. Clustering
        find_clusters(
            obj,
            target='cell',
            k=10,
            res=[0.6, 1.2, 3],
            prefix='pycisTopic_',
            scale=True,
            split_pattern='-'
        )

        # 2. Dimensionality reduction
        run_umap(obj, target='cell', scale=True)
        run_tsne(obj, target='cell', scale=True)

        # 3. Metadata plotting
        variables_to_plot = [cell_type_var] + [f'pycisTopic_leiden_10_{r}' for r in [0.6, 1.2, 3]]
        for var in variables_to_plot:
            plt.figure(figsize=(6,6))
            plot_metadata(
                obj,
                reduction_name='UMAP',
                variables=[var],
                target='cell',
                num_columns=1,
                text_size=10,
                dot_size=5
            )
            output_file = os.path.join(output_dir, f"{sample_name}_{var}_UMAP.pdf")
            plt.savefig(output_file)
            plt.close()
            print(f"Saved plot: {output_file}")

        print(f"Finished {sample_name}\n")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Cluster and visualize pycisTopic MALLET object")
    parser.add_argument("--cistopic_pickle", required=True, help="Pickle file from MALLET with selected models")
    parser.add_argument("--output_dir", required=True, help="Folder to save plots as PDF")
    args = parser.parse_args()

    main(args.cistopic_pickle, args.output_dir)

