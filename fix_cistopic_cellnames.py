#!/usr/bin/env python3
import pickle
import argparse
from collections import defaultdict
import pandas as pd

def clean_cistopic_cellnames(input_pickle, output_pickle):
    """
    Clean CistopicObject cell_names and cell_data.index safely:
    - Remove ___THX suffix
    - Make duplicates unique (_1, _2, ...) if needed
    - Preserve all metadata (celltype_scrna, sample_id, barcode, etc.)
    """
    # Load CistopicObject
    with open(input_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    # --- Step 1: Define cleaning function ---
    def fix_suffix(name):
        """
        Keep everything before ___
        Example:
            GCACTTACAACAGCCT-1-TH1___TH1 -> GCACTTACAACAGCCT-1-TH1
        """
        name = str(name)  # Convert to string to avoid TypeError
        if "___" in name:
            return name.split("___")[0]
        return name

    # --- Step 2: Clean cell_names ---
    if hasattr(cistopic_obj, "cell_names"):
        print("Original first 5 cell_names:", cistopic_obj.cell_names[:5])
        cleaned_cell_names = [fix_suffix(c) for c in cistopic_obj.cell_names]
        cistopic_obj.cell_names = cleaned_cell_names
        print("Cleaned first 5 cell_names:", cistopic_obj.cell_names[:5])

    # --- Step 3: Clean cell_data.index while preserving metadata ---
    if hasattr(cistopic_obj, "cell_data"):
        print("Cleaning cell_data.index to match cleaned cell_names...")

        # Ensure 'barcode' column exists to preserve original cell identifiers
        if 'barcode' not in cistopic_obj.cell_data.columns:
            cistopic_obj.cell_data['barcode'] = cistopic_obj.cell_data.index.astype(str)

        old_index = list(cistopic_obj.cell_data.index)
        cleaned_index = [fix_suffix(c) for c in old_index]

        # Make duplicates unique while preserving all metadata
        counts = defaultdict(int)
        new_index = []
        for c in cleaned_index:
            if counts[c] > 0:
                new_name = f"{c}_{counts[c]}"
            else:
                new_name = c
            counts[c] += 1
            new_index.append(new_name)

        # Assign new index
        cistopic_obj.cell_data = cistopic_obj.cell_data.copy()
        cistopic_obj.cell_data.index = new_index

        # Update cell_names to match cleaned cell_data.index
        cistopic_obj.cell_names = list(cistopic_obj.cell_data.index)

        # Optional: verify metadata is preserved
        if 'celltype_scrna' in cistopic_obj.cell_data.columns:
            print("First 5 celltype_scrna entries:", cistopic_obj.cell_data['celltype_scrna'].head())
        if 'sample_id' in cistopic_obj.cell_data.columns:
            print("First 5 sample_id entries:", cistopic_obj.cell_data['sample_id'].head())
        if 'barcode' in cistopic_obj.cell_data.columns:
            print("First 5 barcode entries:", cistopic_obj.cell_data['barcode'].head())

    # --- Step 4: Save cleaned object ---
    with open(output_pickle, "wb") as f:
        pickle.dump(cistopic_obj, f)

    print(f"Saved cleaned CistopicObject to {output_pickle}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean CistopicObject cell_names and cell_data.index")
    parser.add_argument("--input_pickle", "-i", required=True, help="Input CistopicObject pickle")
    parser.add_argument("--output_pickle", "-o", required=True, help="Output cleaned CistopicObject pickle")
    args = parser.parse_args()

    clean_cistopic_cellnames(args.input_pickle, args.output_pickle)

