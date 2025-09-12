#!/usr/bin/env python3
import os
import argparse
import pandas as pd
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments

def main():
    parser = argparse.ArgumentParser(description="Create cisTopic objects from ATAC fragments.")
    parser.add_argument("--fragments", nargs='+', required=True,
                        help="List of fragment files (TSV.gz) for each sample")
    parser.add_argument("--sample_names", nargs='+', required=True,
                        help="Corresponding sample names")
    parser.add_argument("--consensus_peaks", required=True,
                        help="Consensus peaks BED file")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory to save cistopic objects")
    parser.add_argument("--n_cpu", type=int, default=1,
                        help="Number of CPU cores")
    args = parser.parse_args()

    assert len(args.fragments) == len(args.sample_names), "Number of fragments must match number of sample names"

    os.makedirs(args.output_dir, exist_ok=True)

    for frag_file, sample_name in zip(args.fragments, args.sample_names):
        print(f"Processing sample {sample_name} from {frag_file}...")

        # Prepare metadata: here just a placeholder using barcodes in fragment file
        # For simplicity, pycisTopic can work without extra metadata if you only want peak-by-cell counts
        sample_metadata = pd.DataFrame({"barcode": []})  # empty placeholder, optional

        cistopic_obj = create_cistopic_object_from_fragments(
            path_to_fragments=frag_file,
            path_to_regions=args.consensus_peaks,
            metrics=None,          # optional QC metrics
            valid_bc=None,         # optional valid barcodes
            n_cpu=args.n_cpu,
            project=sample_name,
            split_pattern='-'
        )

        out_file = os.path.join(args.output_dir, f"{sample_name}_cistopic.pkl")
        cistopic_obj.save(out_file)
        print(f"Saved cistopic object for {sample_name} to {out_file}")

if __name__ == "__main__":
    main()

