#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import pyranges as pr
from pycistarget.motif_enrichment_cistarget import cisTarget, cisTargetDatabase

def read_bed_safe(bed_file):
    try:
        if os.path.getsize(bed_file) == 0:
            return None
        df = pd.read_csv(bed_file, sep="\t", header=None, comment='#')
        if df.shape[0] <= 1:
            return None
        if df.shape[1] < 3:
            return None
        df = df.iloc[:, :3]
        df.columns = ['Chromosome', 'Start', 'End']
        return pr.PyRanges(df)
    except:
        return None

def load_region_sets(region_dir):
    region_sets = {}
    for bed_file in os.listdir(region_dir):
        if bed_file.endswith(".bed"):
            path = os.path.join(region_dir, bed_file)
            pr_obj = read_bed_safe(path)
            if pr_obj is not None:
                name = os.path.splitext(bed_file)[0]
                region_sets[name] = pr_obj
    return region_sets

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--feather_file", required=True)
    parser.add_argument("--region_dir", required=True)
    parser.add_argument("--species", default="mus_musculus")
    parser.add_argument("--auc_threshold", type=float, default=0.005)
    parser.add_argument("--nes_threshold", type=float, default=3.0)
    parser.add_argument("--rank_threshold", type=float, default=0.05)
    parser.add_argument("--annotation", nargs="+", default=['Direct_annot', 'Orthology_annot'])
    parser.add_argument("--annotation_version", default="v10nr_clust")
    parser.add_argument("--motif_tbl", required=True)
    parser.add_argument("--output_dir", default="cistarget_results")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    region_sets = load_region_sets(args.region_dir)
    
    if not region_sets:
        print("No valid region sets found")
        return

    # Load database properly
    ctx_db = cisTargetDatabase(args.feather_file)

    results = {}
    for name, region_set in region_sets.items():
        try:
            ct = cisTarget(
                region_set,
                name,
                species=args.species,
                auc_threshold=args.auc_threshold,
                nes_threshold=args.nes_threshold,
                rank_threshold=args.rank_threshold,
                annotation_to_use=args.annotation,
                annotation_version=args.annotation_version,
                path_to_motif_annotations=args.motif_tbl,
            )
            
            result = ct.run_ctx(ctx_db=ctx_db)
            results[name] = result
            
            output_file = os.path.join(args.output_dir, f"{name}_results.csv")
            result.to_csv(output_file)
            print(f"Saved {name} results")
            
        except Exception as e:
            print(f"Error with {name}: {e}")
            continue

    print(f"Completed {len(results)} region sets")

if __name__ == "__main__":
    main()
