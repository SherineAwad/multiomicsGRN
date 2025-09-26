#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import pyranges as pr
from pycistarget.motif_enrichment_cistarget import *

def read_bed_safe(bed_file):
    """Read BED file into PyRanges safely, skipping empty or header-only files."""
    try:
        # First check if file is empty
        if os.path.getsize(bed_file) == 0:
            print(f"Skipping BED {bed_file} (empty file)")
            return None
            
        df = pd.read_csv(bed_file, sep="\t", header=None, comment='#')
        # Skip BEDs with no actual data rows (0 rows)
        if df.shape[0] == 0:
            print(f"Skipping BED {bed_file} (no data rows)")
            return None
            
        # Skip BEDs with only header row (<=1 row after filtering comments)
        if df.shape[0] <= 1:
            print(f"Skipping BED {bed_file} (only header or single row)")
            return None
            
        # Ensure at least 3 columns
        if df.shape[1] < 3:
            print(f"Skipping BED {bed_file} (less than 3 columns)")
            return None
            
        df = df.iloc[:, :3]
        df.columns = ['Chromosome', 'Start', 'End']
        return pr.PyRanges(df)
    except Exception as e:
        print(f"Skipping BED {bed_file} (read error): {e}")
        return None

def load_region_sets(region_dir):
    """Load all BED files in a directory into a dict of PyRanges, skipping empty ones."""
    region_sets = {}
    for bed_file in os.listdir(region_dir):
        if bed_file.endswith(".bed"):
            path = os.path.join(region_dir, bed_file)
            pr_obj = read_bed_safe(path)
            if pr_obj is not None:
                name = os.path.splitext(bed_file)[0]
                region_sets[name] = pr_obj
    print(f"Collected {len(region_sets)} valid region sets.")
    return region_sets

def main():
    parser = argparse.ArgumentParser(description="Run cisTarget on region sets")
    parser.add_argument("--feather_file", required=True, help="Path to cisTarget feather database")
    parser.add_argument("--region_dir", required=True, help="Directory with BED files")
    parser.add_argument("--species", default="mus_musculus", help="Species name")
    parser.add_argument("--auc_threshold", type=float, default=0.005, help="AUC threshold")
    parser.add_argument("--nes_threshold", type=float, default=3.0, help="NES threshold")
    parser.add_argument("--rank_threshold", type=float, default=0.05, help="Rank threshold")
    parser.add_argument("--annotation", nargs="+", default=['Direct_annot', 'Orthology_annot'], help="Annotations")
    parser.add_argument("--annotation_version", default="v10nr_clust", help="Annotation version")
    parser.add_argument("--motif_tbl", required=True, help="Path to motif annotations tbl")
    parser.add_argument("--n_cpu", type=int, default=4, help="Number of CPUs")
    parser.add_argument("--temp_dir", default=None, help="Temporary directory")
    parser.add_argument("--output_dir", default="cistarget_results", help="Output directory for results")
    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Load region sets
    region_sets = load_region_sets(args.region_dir)
    
    # Check if we have any valid region sets
    if len(region_sets) == 0:
        print("ERROR: No valid region sets found. Cannot run cisTarget analysis.")
        return

    # Run cisTarget for each region set separately
    results = {}
    for name, region_set in region_sets.items():
        print(f"Running cisTarget for {name}...")
        
        try:
            # Initialize cistarget object for this region set
            ct = cisTarget(
                region_set,  # First positional argument: the PyRanges object
                name,        # Second positional argument: the name
                species=args.species,  # Corrected from 'specie' to 'species'
                auc_threshold=args.auc_threshold,
                nes_threshold=args.nes_threshold,
                rank_threshold=args.rank_threshold,
                annotation_to_use=args.annotation,  # Corrected parameter name
                annotation_version=args.annotation_version,
                path_to_motif_annotations=args.motif_tbl,
                # Note: n_cpu and temp_dir are not in the signature, so they're removed
            )

            # Run cisTarget
            cistarget_result = ct.run_ctx(ctx_db=args.feather_file)
            results[name] = cistarget_result
            
            # Save results for this region set
            output_file = os.path.join(args.output_dir, f"{name}_cistarget_results.csv")
            # Convert results to DataFrame and save
            if hasattr(cistarget_result, 'to_csv'):
                cistarget_result.to_csv(output_file)
            else:
                # If it's a dictionary or other object, convert to DataFrame first
                pd.DataFrame(cistarget_result).to_csv(output_file)
                
            print(f"Completed {name}. Results saved to {output_file}")
            
        except Exception as e:
            print(f"Error running cisTarget for {name}: {e}")
            continue

    print(f"cisTarget analysis completed. Processed {len(results)} region sets.")
    print(f"Results saved to {args.output_dir}")

if __name__ == "__main__":
    main()
