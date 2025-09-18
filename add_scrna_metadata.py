#!/usr/bin/env python3
import os
import argparse
import pickle
import pandas as pd

def add_scrna_metadata(cistopic_pickle, scrna_csv, output_pickle):
    """
    Add scRNA metadata to a merged CistopicObject.

    cistopic_pickle: pickle containing a merged CistopicObject
    scrna_csv: CSV file with scRNA metadata (columns: 'barcode', 'Seurat_cell_type', etc.)
    output_pickle: path to save updated CistopicObject
    """
    # Load merged CistopicObject
    with open(cistopic_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    # Load scRNA metadata
    scrna_meta = pd.read_csv(scrna_csv)
    scrna_meta.set_index('barcode', inplace=True)

    # Only keep metadata rows for barcodes present in the merged object
    obj_barcodes = cistopic_obj.cell_data.index
    matched_meta = scrna_meta.loc[scrna_meta.index.intersection(obj_barcodes)]

    # Add/merge metadata into CistopicObject.cell_data
    cistopic_obj.cell_data = cistopic_obj.cell_data.join(matched_meta, how='left')

    print(f"Added {matched_meta.shape[1]} scRNA metadata columns for {matched_meta.shape[0]} cells.")

    # Save updated object
    with open(output_pickle, "wb") as f:
        pickle.dump(cistopic_obj, f)

    print(f"Updated merged CistopicObject saved to {output_pickle}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add scRNA metadata to merged CistopicObject")
    parser.add_argument("--cistopic_pickle", required=True,
                        help="Pickle file with merged CistopicObject")
    parser.add_argument("--scrna_csv", required=True,
                        help="CSV file with scRNA metadata (must have 'barcode' column)")
    parser.add_argument("--output_pickle", required=True,
                        help="Output pickle file for updated CistopicObject")
    args = parser.parse_args()

    add_scrna_metadata(args.cistopic_pickle, args.scrna_csv, args.output_pickle)

