#!/usr/bin/env python3
import argparse
import pickle
import numpy as np
import pandas as pd

def filter_cistopic_cells(input_pickle, output_pickle, metadata_col="celltype_scrna"):
    # Load CistopicObject
    with open(input_pickle, "rb") as f:
        obj = pickle.load(f)

    # Check metadata column exists
    if metadata_col not in obj.cell_data.columns:
        raise ValueError(f"Metadata column '{metadata_col}' not found in cell_data")

    # Identify valid cells (non-NaN)
    mask = obj.cell_data[metadata_col].notna()
    num_removed = (~mask).sum()
    print(f"Removing {num_removed} cells with missing '{metadata_col}'")

    # Filter cell_data
    obj.cell_data = obj.cell_data.loc[mask]

    # Filter cell_names to match filtered cell_data
    filtered_cell_names = list(obj.cell_data.index)
    original_cell_names = obj.cell_names.copy()
    obj.cell_names = filtered_cell_names
    print(f"Number of cells remaining: {len(filtered_cell_names)}")

    # Subset cell_topic using integer indices
    if hasattr(obj.selected_model, "cell_topic"):
        cell_topic = obj.selected_model.cell_topic
        if isinstance(cell_topic, pd.DataFrame):
            cell_topic = cell_topic.values  # convert to numpy

        # Map filtered names to column indices in original cell_names
        name_to_idx = {name: i for i, name in enumerate(original_cell_names)}
        col_indices = [name_to_idx[name] for name in filtered_cell_names]

        print(f"Original cell_topic shape: {cell_topic.shape}")
        obj.selected_model.cell_topic = cell_topic[:, col_indices]
        print(f"Filtered cell_topic shape: {obj.selected_model.cell_topic.shape}")

    # Save filtered object
    with open(output_pickle, "wb") as f:
        pickle.dump(obj, f)
    print(f"Filtered object saved to: {output_pickle}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter CistopicObject cells")
    parser.add_argument("-i", "--input_pickle", required=True, help="Input CistopicObject pickle")
    parser.add_argument("-o", "--output_pickle", required=True, help="Output filtered CistopicObject pickle")
    parser.add_argument("-m", "--metadata_col", default="celltype_scrna", help="Metadata column to filter on")
    args = parser.parse_args()
    filter_cistopic_cells(args.input_pickle, args.output_pickle, args.metadata_col)

