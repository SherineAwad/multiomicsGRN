#!/usr/bin/env python3
import os
import argparse
import pickle
from pycisTopic.qc import get_barcodes_passing_qc_for_sample

def main(fragments_dict_path, qc_output_dir, output_pickle):
    """
    Collect barcodes passing QC for each sample.
    
    fragments_dict_path: pickle file containing sample_id -> fragment path
    qc_output_dir: directory where pycisTopic QC outputs are located
    output_pickle: path to save results (barcodes + thresholds)
    """

    # Load fragments_dict
    if not os.path.exists(fragments_dict_path):
        raise FileNotFoundError(f"Fragments dictionary not found: {fragments_dict_path}")

    with open(fragments_dict_path, "rb") as f:
        fragments_dict = pickle.load(f)

    sample_id_to_barcodes_passing_filters = {}
    sample_id_to_thresholds = {}

    # Loop over each sample
    for sample_id in fragments_dict:
        print(f"Processing sample: {sample_id}")
        barcodes, thresholds = get_barcodes_passing_qc_for_sample(
            sample_id=sample_id,
            pycistopic_qc_output_dir=qc_output_dir,
            unique_fragments_threshold=None,  # automatic
            tss_enrichment_threshold=None,    # automatic
            frip_threshold=0,
            use_automatic_thresholds=True
        )
        sample_id_to_barcodes_passing_filters[sample_id] = barcodes
        sample_id_to_thresholds[sample_id] = thresholds

    # Save results
    results = {
        "barcodes": sample_id_to_barcodes_passing_filters,
        "thresholds": sample_id_to_thresholds
    }
    with open(output_pickle, "wb") as f:
        pickle.dump(results, f)

    print(f"QC barcodes and thresholds saved to {output_pickle}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collect barcodes passing pycisTopic QC")
    parser.add_argument("--fragments_dict", required=True,
                        help="Pickle file containing sample_id -> fragment file path")
    parser.add_argument("--qc_output_dir", required=True,
                        help="Directory containing pycisTopic QC output for all samples")
    parser.add_argument("--output_pickle", required=True,
                        help="Path to save barcodes and thresholds pickle")
    args = parser.parse_args()

    main(args.fragments_dict, args.qc_output_dir, args.output_pickle)

