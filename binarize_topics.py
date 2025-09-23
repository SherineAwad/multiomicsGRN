#!/usr/bin/env python3
import argparse
import os
import pickle
from pycisTopic.topic_binarization import binarize_topics
import matplotlib.pyplot as plt
import pandas as pd


def run_binarizations(input_pickle, output_dir):
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
            obj.selected_model.cell_topic.index = [f"Topic{i+1}" for i in range(obj.selected_model.cell_topic.shape[0])]
        else:
            # numpy array
            if cell_topic.shape[1] == len(obj.cell_names):
                obj.selected_model.cell_topic = cell_topic
            obj.selected_model.topic_names = [f"Topic{i+1}" for i in range(cell_topic.shape[0])]

    print("[INFO] LDA model already has cell_topic. Proceeding to binarization...")

    # --- 1) Binarize top n topics ---
    print("[INFO] Running region binarization (top 3000)...")
    region_bin_topics_top_3k = binarize_topics(
        obj, method='ntop', ntop=3000,
        plot=True, num_columns=5
    )
    plt.savefig(os.path.join(output_dir, "region_bin_top3k.png"), bbox_inches='tight')
    plt.close()

    # --- 2) Binarize using Otsu threshold ---
    print("[INFO] Running region binarization (Otsu)...")
    region_bin_topics_otsu = binarize_topics(
        obj, method='otsu',
        plot=True, num_columns=5
    )
    plt.savefig(os.path.join(output_dir, "region_bin_otsu.png"), bbox_inches='tight')
    plt.close()

    # --- 3) Binarize cell topics using Li method ---
    print("[INFO] Running cell binarization (Li)...")
    binarized_cell_topic = binarize_topics(
        obj,
        target='cell',
        method='li',
        plot=True,
        num_columns=5,
        nbins=100
    )
    plt.savefig(os.path.join(output_dir, "cell_topic_li.png"), bbox_inches='tight')
    plt.close()

    # --- 4) Save binarized object ---
    binarized_file = os.path.join(output_dir, "cistopic_obj_binarized.pkl")
    with open(binarized_file, "wb") as f:
        pickle.dump(obj, f)
    print(f"[INFO] Saved binarized cistopic object -> {binarized_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run binarization of CistopicObject topics")
    parser.add_argument("--input_pickle", "-i", required=True, help="Input CistopicObject pickle")
    parser.add_argument("--output_dir", "-o", required=True, help="Directory to save outputs")
    args = parser.parse_args()

    run_binarizations(args.input_pickle, args.output_dir)

