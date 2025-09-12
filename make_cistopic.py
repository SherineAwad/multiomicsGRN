#!/usr/bin/env python3
import os
import argparse
import pandas as pd
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments

def main():
    parser = argparse.ArgumentParser(description="Create cisTopic objects from ATAC fragments and consensus peaks")
    parser.add_argument("--fragments", nargs='+', required=True,
                        help="Paths to fragment files (e.g. TH1_atac_fragments.tsv.gz TH2_atac_fragments.tsv.gz)")
    parser.add_argument("--fragment_names", nargs='+', required=True,
                        help="Names of the corresponding samples (e.g. Control KO)")
    parser.add_argument("--barcodes", required=True,
                        help="CSV file with columns: barcode,celltype,sample")
    parser.add_argument("--consensus_peaks", required=True,
                        help="BED file of consensus peaks")
    parser.add_argument("--blacklist", required=False, default=None,
                        help="BED file of blacklist regions (optional)")
    parser.add_argument("--n_cpu", type=int, default=1,
                        help="Number of CPU threads to use")
    parser.add_argument("--split_pattern", default='-',
                        help="Pattern to split barcodes (default '-')")
    parser.add_argument("--outdir", required=True,
                        help="Directory to save cisTopic objects")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load cell metadata
    cell_data = pd.read_csv(args.barcodes)
    # Strip quotes or spaces from column names
    cell_data.columns = [c.strip().replace('"','') for c in cell_data.columns]

    cistopic_obj_list = []

    if len(args.fragments) != len(args.fragment_names):
        raise ValueError("Number of fragments and fragment_names must match!")

    for frag_file, sample_name in zip(args.fragments, args.fragment_names):
        # Filter metadata for this sample
        sample_metadata = cell_data[cell_data["sample"] == sample_name]
        valid_barcodes = sample_metadata["barcode"].tolist()

        print(f"Creating cisTopic object for sample {sample_name} with {len(valid_barcodes)} cells")

        cistopic_obj = create_cistopic_object_from_fragments(
            path_to_fragments=frag_file,
            path_to_regions=args.consensus_peaks,
            path_to_blacklist=args.blacklist,
            metrics=None,  # Optional QC metrics
            valid_bc=valid_barcodes,
            n_cpu=args.n_cpu,
            project=sample_name,
            split_pattern=args.split_pattern
        )

        # Save cisTopic object as pickle
        out_file = os.path.join(args.outdir, f"cisTopic_{sample_name}.pkl")
        cistopic_obj.save(out_file)
        print(f"Saved cisTopic object to {out_file}")

        cistopic_obj_list.append(cistopic_obj)

    print("All cisTopic objects created successfully.")

if __name__ == "__main__":
    main()

