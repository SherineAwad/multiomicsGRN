#!/usr/bin/env python3
import os
import argparse
import pickle
import numpy as np
import pandas as pd

def filter_empty_scrna_cells(cistopic_pickle, output_pickle, metadata_col='celltype_scrna'):
    """
    Filter out cells with empty or missing scRNA metadata from a CistopicObject
    and save the filtered object.
    """
    # Load CistopicObject(s)
    with open(cistopic_pickle, "rb") as f:
        cistopic_objs = pickle.load(f)

    if not isinstance(cistopic_objs, list):
        cistopic_objs = [cistopic_objs]

    filtered_objs = []
    for idx, obj in enumerate(cistopic_objs):
        print(f"Processing sample {idx}...")
        # Ensure the metadata column exists
        if metadata_col not in obj.cell_data.columns:
            raise ValueError(f"Metadata column '{metadata_col}' not found in CistopicObject")

        # Replace empty strings with NaN
        obj.cell_data[metadata_col] = obj.cell_data[metadata_col].astype(str).replace(r'^\s*$', np.nan, regex=True)

        # Keep only cells with non-empty annotation
        valid_cells = obj.cell_data[metadata_col].notna()
        num_removed = (~valid_cells).sum()
        print(f"Removing {num_removed} cells with empty {metadata_col} out of {obj.cell_data.shape[0]} total")

        # Subset the object
        obj_filtered = obj
        obj_filtered.cell_data = obj_filtered.cell_data.loc[valid_cells]

        print(f"Filtered object has {obj_filtered.cell_data.shape[0]} cells remaining")
        filtered_objs.append(obj_filtered)

    # Save filtered object
    with open(output_pickle, "wb") as f:
        pickle.dump(filtered_objs if len(filtered_objs) > 1 else filtered_objs[0], f)

    print(f"Filtered CistopicObject saved to: {output_pickle}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter empty scRNA cells from CistopicObject")
    parser.add_argument("--cistopic_pickle", required=True, help="Input pickle with merged CistopicObject(s)")
    parser.add_argument("--output_pickle", required=True, help="Output pickle for filtered object")
    parser.add_argument("--metadata_col", default="celltype_scrna", help="Metadata column to filter on (default: celltype_scrna)")
    args = parser.parse_args()

    filter_empty_scrna_cells(args.cistopic_pickle, args.output_pickle, args.metadata_col)

