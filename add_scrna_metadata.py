#!/usr/bin/env python3
import os
import argparse
import pickle
import pandas as pd

def add_scrna_metadata(cistopic_pickle, scrna_csv, output_pickle):
    """
    Add scRNA metadata to a merged CistopicObject.
    - Adds 'celltype' column from scRNA data
    - Removes empty 'celltype_atac' column if it exists
    - Removes cells with missing 'celltype' after join
    """
    # Load CistopicObject
    with open(cistopic_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    # Ensure it's a single object
    if isinstance(cistopic_obj, list):
        if len(cistopic_obj) != 1:
            raise ValueError("Expected a single merged CistopicObject in pickle")
        cistopic_obj = cistopic_obj[0]

    # Load scRNA metadata
    scrna_meta = pd.read_csv(scrna_csv)
    if 'barcode' not in scrna_meta.columns or 'celltype' not in scrna_meta.columns:
        raise ValueError("scRNA CSV must contain 'barcode' and 'celltype' columns")

    # Clean barcodes
    scrna_meta['barcode'] = scrna_meta['barcode'].astype(str).str.strip()
    scrna_meta.set_index('barcode', inplace=True)

    # Prepare ATAC barcodes - barcodes are now in INDEX, not column
    if 'barcode_stripped' not in cistopic_obj.cell_data.columns:
        cistopic_obj.cell_data['barcode_stripped'] = cistopic_obj.cell_data.index.astype(str).str.strip()

    # DEBUG: Check barcode alignment
    print("=== BARCODE DEBUG INFO ===")
    print(f"ATAC barcodes (first 5): {cistopic_obj.cell_data['barcode_stripped'].head(5).tolist()}")
    print(f"scRNA barcodes (first 5): {scrna_meta.index[:5].tolist()}")
    print(f"ATAC shape: {cistopic_obj.cell_data.shape}")
    print(f"scRNA shape: {scrna_meta.shape}")

    cistopic_obj.cell_data.set_index('barcode_stripped', inplace=True)

    # Add scRNA metadata
    sc_col_name = 'celltype'
    if sc_col_name in cistopic_obj.cell_data.columns:
        raise ValueError(f"Column {sc_col_name} already exists")
    matched_meta = scrna_meta[['celltype']].rename(columns={'celltype': sc_col_name})
    cistopic_obj.cell_data = cistopic_obj.cell_data.join(matched_meta, how='left')

    # Report matches
    n_matched = cistopic_obj.cell_data[sc_col_name].notna().sum()
    total_cells = cistopic_obj.cell_data.shape[0]
    print(f"Matched {n_matched} out of {total_cells} cells ({n_matched/total_cells*100:.2f}%)")

    # Remove cells with missing 'celltype'
    before_drop = cistopic_obj.cell_data.shape[0]
    cistopic_obj.cell_data.dropna(subset=[sc_col_name], inplace=True)
    after_drop = cistopic_obj.cell_data.shape[0]
    print(f"✓ Removed {before_drop - after_drop} cells with missing '{sc_col_name}' metadata")

    # Remove empty celltype_atac column if it exists and is all NA
    if 'celltype_atac' in cistopic_obj.cell_data.columns:
        if cistopic_obj.cell_data['celltype_atac'].isna().all():
            cistopic_obj.cell_data.drop(columns=['celltype_atac'], inplace=True)
            print("✓ Removed empty 'celltype_atac' column")
        else:
            print("✓ Kept 'celltype_atac' column (contains some data)")

    # Reset index
    cistopic_obj.cell_data.reset_index(drop=True, inplace=True)

    # Save updated object
    with open(output_pickle, "wb") as f:
        pickle.dump(cistopic_obj, f)
    print(f"✓ Updated CistopicObject saved to: {output_pickle}")
    print(f"✓ Final columns: {list(cistopic_obj.cell_data.columns)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add scRNA metadata to CistopicObject")
    parser.add_argument("--cistopic_pickle", required=True, help="Pickle file with merged CistopicObject")
    parser.add_argument("--scrna_csv", required=True, help="CSV with scRNA metadata (barcode, celltype)")
    parser.add_argument("--output_pickle", required=True, help="Path to save updated pickle")
    args = parser.parse_args()

    add_scrna_metadata(args.cistopic_pickle, args.scrna_csv, args.output_pickle)

