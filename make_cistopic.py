#!/usr/bin/env python3
import os
import argparse
import pandas as pd
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import pickle

def main():
    parser = argparse.ArgumentParser(description="Create cisTopic objects from fragments")
    parser.add_argument("--fragments", nargs='+', required=True,
                        help="Fragment files for each sample (space-separated)")
    parser.add_argument("--barcodes", required=True,
                        help="CSV file with columns: barcode,celltype,sample")
    parser.add_argument("--consensus_bed", required=True,
                        help="Consensus peaks BED file")
    parser.add_argument("--outdir", required=True,
                        help="Output directory to save cisTopic objects")
    parser.add_argument("--cpu", type=int, default=1,
                        help="Number of CPUs to use")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Read barcodes CSV
    cell_data = pd.read_csv(args.barcodes)
    if "barcode" not in cell_data.columns or "sample" not in cell_data.columns:
        raise ValueError("Barcode CSV must have columns: 'barcode', 'celltype', 'sample'")

    # Loop over fragment files and generate cisTopic object
    for frag_path in args.fragments:
        sample_name = os.path.basename(frag_path).split("_")[0]  # e.g., TH1_fragments.tsv.gz -> TH1
        sample_metadata = cell_data[cell_data["sample"] == sample_name]

        print(f"Processing sample: {sample_name}, {len(sample_metadata)} cells")

        cistopic_obj = create_cistopic_object_from_fragments(
            path_to_fragments=frag_path,
            path_to_regions=args.consensus_bed,
            valid_bc=sample_metadata["barcode"].tolist(),
            n_cpu=args.cpu,
            project=sample_name,
            split_pattern="-"
        )

        # Save cisTopic object as pickle
        out_file = os.path.join(args.outdir, f"{sample_name}_cistopic.pkl")
        with open(out_file, "wb") as f:
            pickle.dump(cistopic_obj, f)
        print(f"Saved cisTopic object to: {out_file}")

if __name__ == "__main__":
    main()

