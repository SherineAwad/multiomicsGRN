#!/usr/bin/env python3
import argparse
import pickle
import os
import pandas as pd
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)
import matplotlib.pyplot as plt

def cluster_cistopic(input_pickle, output_pickle, outdir, resolutions=[0.6, 1.2, 3], k=10):
    os.makedirs(outdir, exist_ok=True)

    # Load object
    with open(input_pickle, "rb") as f:
        obj = pickle.load(f)

    print("Original cell_topic shape:", obj.selected_model.cell_topic.shape)
    print("Number of cell_names:", len(obj.cell_names))

    # Ensure cell_topic is a DataFrame with columns matching cell_names
    if not isinstance(obj.selected_model.cell_topic, pd.DataFrame):
        obj.selected_model.cell_topic = pd.DataFrame(
            obj.selected_model.cell_topic,
            columns=obj.cell_names
        )

    print("cell_topic converted to DataFrame:", obj.selected_model.cell_topic.shape)

    # --- Clustering ---
    find_clusters(
        obj,
        target='cell',
        k=k,
        res=resolutions,
        prefix='pycisTopic_',
        scale=True,
        split_pattern='-'
    )

    # --- Dimensionality Reduction ---
    run_umap(obj, target='cell', scale=True)
    run_tsne(obj, target='cell', scale=True)

    # --- Metadata Plots ---
    # Update this list according to your metadata columns
    meta_vars = ['celltype_scrna'] + [f'pycisTopic_leiden_{k}_{r}' for r in resolutions]
    plot_metadata(
        obj,
        reduction_name='UMAP',
        variables=meta_vars,
        target='cell',
        num_columns=4,
        text_size=10,
        dot_size=5
    )
    plt.savefig(os.path.join(outdir, "metadata_umap.png"), bbox_inches='tight')
    plt.close()

    # --- Annotate clusters with most frequent cell type ---
    annot_dict = {}
    for res in resolutions:
        col = f'pycisTopic_leiden_{k}_{res}'
        annot_dict[col] = {}
        clusters = set(obj.cell_data[col])
        for cluster in clusters:
            idx = obj.cell_data[obj.cell_data[col] == cluster].index
            counts = obj.cell_data.loc[idx, 'celltype_scrna'].value_counts()
            annot_dict[col][cluster] = f"{counts.index[0]}({cluster})"
        # Apply annotations
        obj.cell_data[col] = [annot_dict[col][x] for x in obj.cell_data[col].tolist()]

    # Plot annotated clusters
    plot_metadata(
        obj,
        reduction_name='UMAP',
        variables=[f'pycisTopic_leiden_{k}_{r}' for r in resolutions],
        target='cell',
        num_columns=3,
        text_size=10,
        dot_size=5
    )
    plt.savefig(os.path.join(outdir, "annotated_clusters_umap.png"), bbox_inches='tight')
    plt.close()

    # --- Additional QC plots if present ---
    qc_vars = ['log10_unique_fragments_count', 'tss_enrichment', 'Doublet_scores_fragments', 'fraction_of_fragments_in_peaks']
    existing_qc_vars = [v for v in qc_vars if v in obj.cell_data.columns]
    if existing_qc_vars:
        plot_metadata(
            obj,
            reduction_name='UMAP',
            variables=existing_qc_vars,
            target='cell',
            num_columns=4,
            text_size=10,
            dot_size=5
        )
        plt.savefig(os.path.join(outdir, "qc_metrics_umap.png"), bbox_inches='tight')
        plt.close()

    # --- Topic Plots ---
    plot_topic(
        obj,
        reduction_name='UMAP',
        target='cell',
        num_columns=5
    )
    plt.savefig(os.path.join(outdir, "topic_umap.png"), bbox_inches='tight')
    plt.close()

    # --- Cell-topic heatmap ---
    cell_topic_heatmap(
        obj,
        variables=['celltype_scrna'],
        scale=False,
        legend_loc_x=1.0,
        legend_loc_y=-1.2,
        legend_dist_y=-1,
        figsize=(10, 10)
    )
    plt.savefig(os.path.join(outdir, "cell_topic_heatmap.png"), bbox_inches='tight')
    plt.close()

    # Save filtered and clustered object
    with open(output_pickle, "wb") as f:
        pickle.dump(obj, f)

    print(f"Filtered and clustered CistopicObject saved to: {output_pickle}")
    print("All plots saved to:", outdir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster, UMAP, annotate, and plot a CistopicObject")
    parser.add_argument("-i", "--input_pickle", required=True, help="Input CistopicObject pickle")
    parser.add_argument("-o", "--output_pickle", required=True, help="Output clustered CistopicObject pickle")
    parser.add_argument("-d", "--outdir", required=True, help="Directory to save plots")
    parser.add_argument("--resolutions", nargs="+", type=float, default=[0.6, 1.2, 3], help="Leiden resolutions")
    parser.add_argument("--k", type=int, default=10, help="Number of neighbors for clustering")

    args = parser.parse_args()
    cluster_cistopic(args.input_pickle, args.output_pickle, args.outdir, args.resolutions, args.k)

