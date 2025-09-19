#!/usr/bin/env python3
import pickle
import argparse
import numpy as np
import pandas as pd
import os

def explore_cistopic_pickle(pkl_file):
    with open(pkl_file, "rb") as f:
        obj = pickle.load(f)

    # Make it always a list
    if not isinstance(obj, list):
        obj = [obj]

    print(f"Total CistopicObjects in pickle: {len(obj)}\n")

    for idx, cistopic_obj in enumerate(obj):
        print(f"--- Sample {idx} ---")

        # 1. Cell counts
        n_cells = cistopic_obj.cell_data.shape[0] if hasattr(cistopic_obj, "cell_data") else "N/A"
        print(f"Number of cells: {n_cells}")

        # 2. Peak counts
        n_peaks = cistopic_obj.counts.shape[1] if hasattr(cistopic_obj, "counts") else "N/A"
        print(f"Number of peaks: {n_peaks}")

        # 3. Fragment distribution
        if hasattr(cistopic_obj, "cell_data") and "total_fragments_count" in cistopic_obj.cell_data.columns:
            counts = cistopic_obj.cell_data["total_fragments_count"]
            print("Fragment counts per cell (summary):")
            print(f"  min: {counts.min()}, max: {counts.max()}, median: {counts.median()}, mean: {counts.mean():.1f}")
            low_cells = (counts < 1000).sum()
            print(f"  Cells with <1000 fragments: {low_cells}")
        else:
            print("No fragment count info available")

        # 4. Check selected_model
        if hasattr(cistopic_obj, "selected_model"):
            sel_model = cistopic_obj.selected_model
            if sel_model is None:
                print("selected_model: None")
            elif isinstance(sel_model, list):
                print(f"selected_model: list of {len(sel_model)} models")
                for m_idx, m in enumerate(sel_model):
                    n_topics = len(m.topic_names) if hasattr(m, "topic_names") else "N/A"
                    print(f"  Model {m_idx}: {n_topics} topics")
                    if hasattr(m, "cell_topic") and m.cell_topic is not None:
                        print(f"    cell_topic shape: {m.cell_topic.shape}")
                    else:
                        print("    cell_topic: None")
            else:
                n_topics = len(sel_model.topic_names) if hasattr(sel_model, "topic_names") else "N/A"
                print(f"selected_model: single model with {n_topics} topics")
                if hasattr(sel_model, "cell_topic") and sel_model.cell_topic is not None:
                    print(f"  cell_topic shape: {sel_model.cell_topic.shape}")
                else:
                    print("  cell_topic: None")
        else:
            print("No selected_model attribute")

        # 5. Metadata keys
        if hasattr(cistopic_obj, "cell_data"):
            print(f"Metadata keys: {list(cistopic_obj.cell_data.keys())}\n")
        else:
            print("No cell_data attribute\n")

        # 6. Check for clustering columns starting with pycisTopic_leiden
        if hasattr(cistopic_obj, "cell_data"):
            cluster_cols = [col for col in cistopic_obj.cell_data.columns if col.startswith("pycisTopic_leiden")]
            if cluster_cols:
                print("Found clustering columns:")
                for col in cluster_cols:
                    n_clusters = cistopic_obj.cell_data[col].nunique()
                    print(f"  {col}: {n_clusters} clusters")
            else:
                print("No pycisTopic_leiden clustering columns found.")

        # 7. Check counts attribute
        print("Has counts:", hasattr(cistopic_obj, "counts"))
        print("Counts type:", type(getattr(cistopic_obj, "counts", None)))
        print("\n" + "-"*50 + "\n")

        # 8. Quick summaries
        if hasattr(cistopic_obj, "cell_data"):
            if 'sample_id' in cistopic_obj.cell_data.columns:
                print("sample id", cistopic_obj.cell_data['sample_id'].value_counts())
            if 'celltype_atac' in cistopic_obj.cell_data.columns:
                print("Celltype ATAC", cistopic_obj.cell_data['celltype_atac'].value_counts())
            if 'celltype_scrna' in cistopic_obj.cell_data.columns:
                print("Celltype scRNA (top 10):")
                print(cistopic_obj.cell_data['celltype_scrna'].value_counts().head(10))
                print(cistopic_obj.cell_data['celltype_scrna'].unique())
                print(cistopic_obj.cell_data['celltype_scrna'].iloc[:20])
                print("Number of NaN in celltype_scrna:", cistopic_obj.cell_data['celltype_scrna'].isna().sum())

                # Save all scRNA celltypes to a file
                output_file = "scRNA_cistopic.out"
                cistopic_obj.cell_data['celltype_scrna'].to_csv(
                    output_file,
                    index=True,
                    header=True,
                    sep="\t"
                )
                print(f"Saved full celltype_scrna column to {output_file}")

        # 9. Transpose cell_topic (optional, keep your code)
        if hasattr(cistopic_obj, "selected_model") and hasattr(cistopic_obj.selected_model, "cell_topic"):
            meta = cistopic_obj.cell_data
            cell_topic = pd.DataFrame(
                cistopic_obj.selected_model.cell_topic.T,
                index=cistopic_obj.cell_names,
                columns=[f"Topic_{i}" for i in range(cistopic_obj.selected_model.cell_topic.shape[0])]
            )
            # Fix KeyError by selecting intersection of indices
            common_index = cell_topic.index.intersection(meta.index)
            cell_topic = cell_topic.loc[common_index]

        # --- APPENDED SIMPLE ATTRIBUTE EXPLORATION ---
        print("\nTop-level attributes of the CistopicObject:")
        for attr in dir(cistopic_obj):
            if not attr.startswith("__"):
                try:
                    val = getattr(cistopic_obj, attr)
                    summary = ""
                    if isinstance(val, (list, tuple, dict, set)):
                        summary = f"len={len(val)}"
                    elif hasattr(val, "shape"):
                        summary = f"shape={val.shape}"
                    elif isinstance(val, str):
                        summary = f"value='{val[:50]}'"
                    print(f"  {attr} ({type(val)}) {summary}")
                except Exception:
                    print(f"  {attr} ({type(val)}) - could not inspect")
        print("\n" + "="*80 + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Explore CistopicObject pickle")
    parser.add_argument("--input", "-i", required=True, help="Input pickle file")
    args = parser.parse_args()

    explore_cistopic_pickle(args.input)

