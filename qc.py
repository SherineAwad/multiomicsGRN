#!/usr/bin/env python3
import os
import gzip
import argparse
import pickle
from collections import defaultdict

def read_fragments(fragment_file):
    """Return barcode -> fragment counts dict"""
    barcode_counts = defaultdict(int)
    with gzip.open(fragment_file, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            barcode = parts[3]
            barcode_counts[barcode] += 1
    return barcode_counts

def apply_threshold(barcode_counts, threshold):
    """Return barcodes passing unique fragment threshold"""
    return [b for b, c in barcode_counts.items() if c >= threshold]

def main():
    parser = argparse.ArgumentParser(description="QC for multiple ATAC-seq fragment files")
    parser.add_argument("--fragments", nargs='+', required=True, help="Fragment files (*.tsv.gz)")
    parser.add_argument("--unique_fragments_threshold", type=int, default=None, help="Custom unique fragment cutoff")
    parser.add_argument("--output_dir", required=True, help="Directory to save QC pickles")
    parser.add_argument("--n_top", type=int, default=5, help="Top N barcodes to show per sample")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    qc_default_dir = os.path.join(args.output_dir, "QC_default")
    qc_custom_dir = os.path.join(args.output_dir, "QC_custom")
    os.makedirs(qc_default_dir, exist_ok=True)
    os.makedirs(qc_custom_dir, exist_ok=True)

    default_qc = {'barcodes': {}, 'thresholds': {}}
    custom_qc = {'barcodes': {}, 'thresholds': {}}

    for fpath in args.fragments:
        sample = os.path.basename(fpath).replace("_atac_fragments.tsv.gz","")
        print(f"\n=== Processing {fpath} ({sample}) ===")

        barcode_counts = read_fragments(fpath)
        total_barcodes = len(barcode_counts)
        avg_fragments = sum(barcode_counts.values()) / total_barcodes if total_barcodes else 0

        # --- Default QC (PyCistopic-style, threshold = average fragments) ---
        default_threshold = int(avg_fragments)
        kept_default = apply_threshold(barcode_counts, default_threshold)
        default_qc['barcodes'][sample] = kept_default
        default_qc['thresholds'][sample] = {'unique_fragments_threshold': default_threshold}

        print(f"[Default QC] Total barcodes: {total_barcodes}, Kept after QC: {len(kept_default)}, unique_fragments_threshold: {default_threshold}")
        top_default = sorted(barcode_counts.items(), key=lambda x: x[1], reverse=True)[:args.n_top]
        print(f"Top {args.n_top} barcodes after default QC:")
        for b, c in top_default:
            print(f"   {b}: {c} fragments")

        # --- Custom QC if specified ---
        if args.unique_fragments_threshold is not None:
            kept_custom = apply_threshold(barcode_counts, args.unique_fragments_threshold)
            custom_qc['barcodes'][sample] = kept_custom
            custom_qc['thresholds'][sample] = {'unique_fragments_threshold': args.unique_fragments_threshold}

            print(f"[Manual QC] Total barcodes: {total_barcodes}, Kept after QC: {len(kept_custom)}, unique_fragments_threshold: {args.unique_fragments_threshold}")
            top_custom = sorted(barcode_counts.items(), key=lambda x: x[1], reverse=True)[:args.n_top]
            print(f"Top {args.n_top} barcodes after manual QC:")
            for b, c in top_custom:
                print(f"   {b}: {c} fragments")

    # --- Save pickles ---
    default_pickle = os.path.join(qc_default_dir, "qc_barcodes_thresholds.pkl")
    custom_pickle = os.path.join(qc_custom_dir, "qc_barcodes_thresholds.pkl")

    with open(default_pickle, "wb") as f:
        pickle.dump(default_qc, f)
    print(f"\n✅ Default QC pickle saved to: {default_pickle}")

    if args.unique_fragments_threshold is not None:
        with open(custom_pickle, "wb") as f:
            pickle.dump(custom_qc, f)
        print(f"✅ Custom QC pickle saved to: {custom_pickle}")

if __name__ == "__main__":
    main()

