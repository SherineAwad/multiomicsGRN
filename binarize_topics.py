#!/usr/bin/env python3
import argparse
import os
import pickle
import numpy as np
import pandas as pd
from pycisTopic.topic_binarization import binarize_topics
import matplotlib.pyplot as plt


def run_binarizations(input_pickle, output_dir, ntop):
    os.makedirs(output_dir, exist_ok=True)

    # Load CistopicObject
    with open(input_pickle, "rb") as f:
        obj = pickle.load(f)

    # Add selected LDA model if missing
    if hasattr(obj, "selected_model") and obj.selected_model is not None:
        if not hasattr(obj.selected_model, "cell_topic") or obj.selected_model.cell_topic is None:
            print("[INFO] Adding selected_model to object...")
            obj.add_LDA_model(obj.selected_model)
    else:
        raise ValueError("[ERROR] No selected_model found in object. Add LDA model first.")

    # Ensure cell_topic has proper topic names
    if hasattr(obj.selected_model, 'cell_topic'):
        cell_topic = obj.selected_model.cell_topic
        if isinstance(cell_topic, pd.DataFrame):
            # Transpose if index are cells
            if list(cell_topic.index) == obj.cell_names:
                obj.selected_model.cell_topic = cell_topic.T
            # Assign Topic1, Topic2 ...
            obj.selected_model.cell_topic.index = [
                f"Topic{i+1}" for i in range(obj.selected_model.cell_topic.shape[0])
            ]
        else:
            # numpy array
            if cell_topic.shape[1] == len(obj.cell_names):
                obj.selected_model.cell_topic = cell_topic
            obj.selected_model.topic_names = [
                f"Topic{i+1}" for i in range(cell_topic.shape[0])
            ]

    print("[INFO] LDA model already has cell_topic. Proceeding to binarization...")

    # --- 1) Binarize top n topics ---
    print(f"[INFO] Running region binarization (top {ntop})...")
    region_bin_topics_top_dict = binarize_topics(
        obj, method='ntop', ntop=ntop,
        plot=True, num_columns=5
    )
    plt.savefig(os.path.join(output_dir, f"region_bin_top{ntop}.png"), bbox_inches='tight')
    plt.close()

    # Convert dictionary to combined matrix
    region_bin_topics_top_combined = _combine_topic_matrices(region_bin_topics_top_dict, obj.region_names)
    print(f"   Combined matrix shape: {region_bin_topics_top_combined.shape}")
    print(f"   Non-zero entries: {(region_bin_topics_top_combined > 0).sum().sum()}")

    # --- 2) Binarize using Otsu threshold ---
    print("[INFO] Running region binarization (Otsu)...")
    region_bin_topics_otsu_dict = binarize_topics(
        obj, method='otsu',
        plot=True, num_columns=5
    )
    plt.savefig(os.path.join(output_dir, "region_bin_otsu.png"), bbox_inches='tight')
    plt.close()

    # Convert dictionary to combined matrix
    region_bin_topics_otsu_combined = _combine_topic_matrices(region_bin_topics_otsu_dict, obj.region_names)
    print(f"   Combined matrix shape: {region_bin_topics_otsu_combined.shape}")
    print(f"   Non-zero entries: {(region_bin_topics_otsu_combined > 0).sum().sum()}")

    # --- 3) Binarize cell topics using Li method ---
    print("[INFO] Running cell binarization (Li)...")
    binarized_cell_topic_dict = binarize_topics(
        obj,
        target='cell',
        method='li',
        plot=True,
        num_columns=5,
        nbins=100
    )
    plt.savefig(os.path.join(output_dir, "cell_topic_li.png"), bbox_inches='tight')
    plt.close()

    # Convert dictionary to combined matrix
    binarized_cell_topic_combined = _combine_cell_matrices(binarized_cell_topic_dict, obj.cell_names)
    print(f"   Combined matrix shape: {binarized_cell_topic_combined.shape}")
    print(f"   Non-zero entries: {(binarized_cell_topic_combined > 0).sum().sum()}")

    # --- Attach binarized outputs to the object ---
    # Use the COMBINED matrices, not the dictionaries
    obj.topic_binarized_matrix = region_bin_topics_top_combined  # This is what dAR needs
    obj.region_bin_topics_top = region_bin_topics_top_combined
    obj.region_bin_topics_otsu = region_bin_topics_otsu_combined
    obj.binarized_cell_topic = binarized_cell_topic_combined

    # --- 4) Save full binarized object ---
    binarized_file = os.path.join(output_dir, "cistopic_obj_binarized.pkl")
    with open(binarized_file, "wb") as f:
        pickle.dump(obj, f)
    print(f"[INFO] Saved binarized cistopic object -> {binarized_file}")

    # --- 5) Save binarized cell_topic DataFrame separately ---
    binarized_cell_topic_file = os.path.join(output_dir, "binarized_cell_topic.pkl")
    with open(binarized_cell_topic_file, "wb") as f:
        pickle.dump(binarized_cell_topic_combined, f)
    print(f"[INFO] Saved binarized cell_topic DataFrame -> {binarized_cell_topic_file}")


def _combine_topic_matrices(topic_dict, all_region_names):
    """
    Combine individual topic binary matrices into one combined matrix
    """
    n_regions = len(all_region_names)
    n_topics = len(topic_dict)
    topic_names = list(topic_dict.keys())

    # Create empty matrix with all regions
    combined_matrix = np.zeros((n_regions, n_topics), dtype=int)
    combined_df = pd.DataFrame(combined_matrix,
                              index=all_region_names,
                              columns=topic_names)

    # Fill in the selected regions for each topic
    for topic_name, topic_matrix in topic_dict.items():
        if isinstance(topic_matrix, pd.DataFrame):
            # Get regions that are selected (value = 1)
            # Since binarize_topics returns probabilities, we need to check which regions are actually selected
            # For top-n method, all regions in the matrix are the selected ones
            selected_regions = topic_matrix.index
            combined_df.loc[selected_regions, topic_name] = 1

    return combined_df


def _combine_cell_matrices(topic_dict, all_cell_names):
    """
    Combine individual cell binary matrices into one combined matrix
    """
    n_cells = len(all_cell_names)
    n_topics = len(topic_dict)
    topic_names = list(topic_dict.keys())

    # Create empty matrix with all cells
    combined_matrix = np.zeros((n_topics, n_cells), dtype=int)
    combined_df = pd.DataFrame(combined_matrix,
                              index=topic_names,
                              columns=all_cell_names)

    # Fill in the selected cells for each topic
    for topic_name, topic_matrix in topic_dict.items():
        if isinstance(topic_matrix, pd.DataFrame):
            # Get cells that are selected (value = 1)
            # For cell binarization, cells are usually the index
            selected_cells = topic_matrix.index
            combined_df.loc[topic_name, selected_cells] = 1

    return combined_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run binarization of CistopicObject topics")
    parser.add_argument("--input_pickle", "-i", required=True, help="Input CistopicObject pickle")
    parser.add_argument("--output_dir", "-o", required=True, help="Directory to save outputs")
    parser.add_argument("--ntop", "-n", type=int, default=10000, help="Number of top regions to select (default: 10000)")
    args = parser.parse_args()

    run_binarizations(args.input_pickle, args.output_dir, args.ntop)
