#!/usr/bin/env python3
import os
import argparse
import pickle
from pycisTopic.plotting.qc_plot import plot_barcode_stats
import matplotlib.pyplot as plt

def main(fragments_dict_path, qc_results_pickle, qc_output_dir, plots_output_dir):
    """
    Generate barcode-level QC plots for each sample and save as PNG.
    
    fragments_dict_path: pickle file containing sample_id -> fragment path
    qc_results_pickle: pickle file with barcodes + thresholds (from QC filtering)
    qc_output_dir: directory containing pycisTopic QC outputs
    plots_output_dir: directory to save QC plots (PNG)
    """

    # Load fragments_dict
    if not os.path.exists(fragments_dict_path):
        raise FileNotFoundError(f"Fragments dictionary not found: {fragments_dict_path}")
    with open(fragments_dict_path, "rb") as f:
        fragments_dict = pickle.load(f)

    # Load QC results (barcodes + thresholds)
    if not os.path.exists(qc_results_pickle):
        raise FileNotFoundError(f"QC results pickle not found: {qc_results_pickle}")
    with open(qc_results_pickle, "rb") as f:
        qc_results = pickle.load(f)

    sample_id_to_barcodes_passing_filters = qc_results["barcodes"]
    sample_id_to_thresholds = qc_results["thresholds"]

    # Create plots output directory
    os.makedirs(plots_output_dir, exist_ok=True)

    # Generate plots
    for sample_id in fragments_dict:
        print(f"Plotting barcode QC for sample: {sample_id}")
        fig = plot_barcode_stats(
            sample_id=sample_id,
            pycistopic_qc_output_dir=qc_output_dir,
            bc_passing_filters=sample_id_to_barcodes_passing_filters[sample_id],
            detailed_title=False,
            **sample_id_to_thresholds[sample_id]
        )

        # Save figure
        plot_file = os.path.join(plots_output_dir, f"{sample_id}_barcode_qc.png")
        fig.savefig(plot_file, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved barcode QC plot: {plot_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate barcode-level pycisTopic QC plots")
    parser.add_argument("--fragments_dict", required=True,
                        help="Pickle file containing sample_id -> fragment file path")
    parser.add_argument("--qc_results_pickle", required=True,
                        help="Pickle file containing barcodes + thresholds from QC filtering")
    parser.add_argument("--qc_output_dir", required=True,
                        help="Directory containing pycisTopic QC output for all samples")
    parser.add_argument("--plots_output_dir", required=True,
                        help="Directory to save barcode QC plots as PNG")
    args = parser.parse_args()

    main(args.fragments_dict, args.qc_results_pickle, args.qc_output_dir, args.plots_output_dir)

