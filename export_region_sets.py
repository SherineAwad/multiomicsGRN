#!/usr/bin/env python3
import os
import pickle
import argparse
from pycisTopic.utils import region_names_to_coordinates

def export_region_sets(cistopic_pickle, out_dir):
    # Create folders
    os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "region_sets", "DARs_cell_type"), exist_ok=True)

    # Load binarised CistopicObject
    with open(cistopic_pickle, "rb") as f:
        obj = pickle.load(f)

    # --- Export topics region sets ---
    if hasattr(obj, "region_bin_topics_otsu"):
        for topic in obj.region_bin_topics_otsu:
            df = region_names_to_coordinates(obj.region_bin_topics_otsu[topic].index)
            df.sort_values(["Chromosome", "Start", "End"]).to_csv(
                os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
                sep="\t", header=False, index=False
            )
    else:
        print("[WARNING] region_bin_topics_otsu not found in object. Skipping Topics_otsu export.")

    if hasattr(obj, "region_bin_topics_top_3k"):
        for topic in obj.region_bin_topics_top_3k:
            df = region_names_to_coordinates(obj.region_bin_topics_top_3k[topic].index)
            df.sort_values(["Chromosome", "Start", "End"]).to_csv(
                os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
                sep="\t", header=False, index=False
            )
    else:
        print("[WARNING] region_bin_topics_top_3k not found in object. Skipping Topics_top_3k export.")

    # --- Export DARs per cell type ---
    if hasattr(obj, "markers_dict"):
        for cell_type in obj.markers_dict:
            df = region_names_to_coordinates(obj.markers_dict[cell_type].index)
            df.sort_values(["Chromosome", "Start", "End"]).to_csv(
                os.path.join(out_dir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
                sep="\t", header=False, index=False
            )
    else:
        print("[WARNING] markers_dict not found in object. Skipping DARs export.")

    print(f"[INFO] Region sets exported to: {os.path.join(out_dir, 'region_sets')}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export Topics and DARs region sets from a binarised CistopicObject")
    parser.add_argument("-i", "--input_pickle", required=True, help="Binarised CistopicObject pickle file")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for BED files")
    args = parser.parse_args()

    export_region_sets(args.input_pickle, args.out_dir)

