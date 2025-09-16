#!/usr/bin/env python3
import argparse
import pickle
import os
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)
import matplotlib.pyplot as plt

def main(args):
    # Load input cisTopic object
    with open(args.input, "rb") as f:
        cistopic_obj = pickle.load(f)

    # Run clustering
    cistopic_obj = find_clusters(
        cistopic_obj,
        target="cell",
        k=args.k,
        res=args.res,
        prefix=args.prefix,
        scale=True,
        split_pattern="-"
    )

    # Run UMAP only
    cistopic_obj = run_umap(cistopic_obj, target="cell")

    # Create output directory for plots
    os.makedirs(args.plot_dir, exist_ok=True)

    # Plot metadata
    fig = plot_metadata(
        cistopic_obj,
        reduction_name="UMAP",
        variables=args.metadata_vars,
        target="cell",
        num_columns=4,
        text_size=10,
        dot_size=5
    )
    fig.savefig(os.path.join(args.plot_dir, "metadata_umap.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Example: topic plot (first topic only, you can expand as needed)
    fig = plot_topic(cistopic_obj, reduction_name="UMAP", topic=1, target="cell")
    fig.savefig(os.path.join(args.plot_dir, "topic_umap.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Example: cell-topic heatmap
    fig = cell_topic_heatmap(cistopic_obj, scale=True)
    fig.savefig(os.path.join(args.plot_dir, "cell_topic_heatmap.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Save updated cisTopic object
    with open(args.output, "wb") as f:
        pickle.dump(cistopic_obj, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run clustering, UMAP, and plotting with pycisTopic"
    )
    parser.add_argument(
        "--input", "-i", required=True, help="Input cisTopic object (pickle file)"
    )
    parser.add_argument(
        "--output", "-o", required=True, help="Output cisTopic object (pickle file)"
    )
    parser.add_argument(
        "--plot_dir", "-p", required=True, help="Directory to save plots"
    )
    parser.add_argument(
        "--k", type=int, default=10, help="k parameter for nearest neighbors"
    )
    parser.add_argument(
        "--res",
        type=float,
        nargs="+",
        default=[0.6, 1.2, 3.0],
        help="Resolution values for clustering"
    )
    parser.add_argument(
        "--prefix", default="pycisTopic_", help="Prefix for cluster labels"
    )
    parser.add_argument(
        "--metadata_vars",
        nargs="+",
        default=["Seurat_cell_type", "pycisTopic_leiden_10_0.6", "pycisTopic_leiden_10_1.2", "pycisTopic_leiden_10_3"],
        help="Metadata variables to plot on UMAP"
    )

    args = parser.parse_args()
    main(args)

