#!/usr/bin/env python3
import os
import argparse
import pickle
from pycisTopic.plotting.qc_plot import plot_sample_stats
import matplotlib.pyplot as plt

def main(fragments_dict_path, qc_output_dir, plots_output_dir):
    """
    Generate QC plots for each sample and save as PNG.
    
    fragments_dict_path: pickle file containing sample_id -> fragment path
    qc_output_dir: directory where pycisTopic QC outputs are located
    plots_output_dir: directory to save QC plots (PNG)
    """

    # Load fragments_dict
    if not os.path.exists(fragments_dict_path):
        raise FileNotFoundError(f"Fragments dictionary not found: {fragments_dict_path}")

    with open(fragments_dict_path, "rb") as f:
        fragments_dict = pickle.load(f)

    # Create plots output directory
    os.makedirs(plots_output_dir, exist_ok=True)

    for sample_id in fragments_dict:
        print(f"Plotting QC for sample: {sample_id}")
        fig = plot_sample_stats(
            sample_id=sample_id,
            pycistopic_qc_output_dir=qc_output_dir
        )

        # Save figure as PNG
        plot_file = os.path.join(plots_output_dir, f"{sample_id}_qc.png")
        fig.savefig(plot_file, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved plot: {plot_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate and save pycisTopic QC plots")
    parser.add_argument("--fragments_dict", required=True,
                        help="Pickle file containing sample_id -> fragment file path")
    parser.add_argument("--qc_output_dir", required=True,
                        help="Directory containing pycisTopic QC output for all samples")
    parser.add_argument("--plots_output_dir", required=True,
                        help="Directory to save QC plots as PNG")
    args = parser.parse_args()

    main(args.fragments_dict, args.qc_output_dir, args.plots_output_dir)

