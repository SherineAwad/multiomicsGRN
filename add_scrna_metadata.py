#!/usr/bin/env python3
import os
import argparse
import pickle
import pandas as pd

def add_scrna_metadata(cistopic_pickle, scrna_csv, output_pickle):
    """
    Add scRNA metadata to a merged CistopicObject and preserve ATAC celltypes.

    - Existing 'celltype' in ATAC will be renamed to 'celltype_atac'.
    - scRNA metadata will be added as 'celltype_scrna'.
    """
    # Load CistopicObject
    with open(cistopic_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    # Ensure it's a single object
    if isinstance(cistopic_obj, list):
        if len(cistopic_obj) != 1:
            raise ValueError("Expected a single merged CistopicObject in pickle")
        cistopic_obj = cistopic_obj[0]

    # Rename existing celltype if exists
    if 'celltype' in cistopic_obj.cell_data.columns:
        if 'celltype_atac' in cistopic_obj.cell_data.columns:
            raise ValueError("Column 'celltype_atac' already exists")
        cistopic_obj.cell_data.rename(columns={'celltype': 'celltype_atac'}, inplace=True)
        print("Renamed existing 'celltype' to 'celltype_atac'")
    else:
        cistopic_obj.cell_data['celltype_atac'] = pd.NA

    # Load scRNA metadata
    scrna_meta = pd.read_csv(scrna_csv)
    if 'barcode' not in scrna_meta.columns or 'celltype' not in scrna_meta.columns:
        raise ValueError("scRNA CSV must contain 'barcode' and 'celltype' columns")

    # Clean barcodes
    scrna_meta['barcode'] = scrna_meta['barcode'].astype(str).str.strip()
    scrna_meta.set_index('barcode', inplace=True)

    # Prepare ATAC barcodes
    if 'barcode_stripped' not in cistopic_obj.cell_data.columns:
        cistopic_obj.cell_data['barcode_stripped'] = cistopic_obj.cell_data['barcode'].astype(str).str.strip()

    cistopic_obj.cell_data.set_index('barcode_stripped', inplace=True)

    # Add scRNA metadata
    sc_col_name = 'celltype_scrna'
    if sc_col_name in cistopic_obj.cell_data.columns:
        raise ValueError(f"Column {sc_col_name} already exists")
    matched_meta = scrna_meta[['celltype']].rename(columns={'celltype': sc_col_name})
    cistopic_obj.cell_data = cistopic_obj.cell_data.join(matched_meta, how='left')

    # Report matches
    n_matched = cistopic_obj.cell_data[sc_col_name].notna().sum()
    total_cells = cistopic_obj.cell_data.shape[0]
    print(f"Matched {n_matched} out of {total_cells} cells ({n_matched/total_cells*100:.2f}%)")

    # Reset index
    cistopic_obj.cell_data.reset_index(drop=True, inplace=True)

    # Save updated object
    with open(output_pickle, "wb") as f:
        pickle.dump(cistopic_obj, f)
    print(f"Updated CistopicObject saved to: {output_pickle}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add scRNA metadata and preserve ATAC celltypes")
    parser.add_argument("--cistopic_pickle", required=True, help="Pickle file with merged CistopicObject")
    parser.add_argument("--scrna_csv", required=True, help="CSV with scRNA metadata (barcode, celltype)")
    parser.add_argument("--output_pickle", required=True, help="Path to save updated pickle")
    args = parser.parse_args()

    add_scrna_metadata(args.cistopic_pickle, args.scrna_csv, args.output_pickle)

