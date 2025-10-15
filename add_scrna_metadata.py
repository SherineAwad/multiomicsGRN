#!/usr/bin/env python3
import os
import argparse
import pickle
import pandas as pd
import re

def add_scrna_metadata(cistopic_pickle, scrna_csv, output_pickle):
    """
    Add scRNA metadata (celltype) to a merged cisTopicObject.
    Includes:
      - Barcode cleaning and alignment
      - Metadata join
      - Removal of unmatched cells
      - Synchronization of .cell_names with filtered cells
      - Optional cleanup of weird suffixes (-TH1___TH1 etc.)
    """
    print(f"🔹 Loading cisTopic object: {cistopic_pickle}")
    with open(cistopic_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    # --- Load and clean scRNA metadata ---
    scrna_meta = pd.read_csv(scrna_csv)
    if 'barcode' not in scrna_meta.columns or 'celltype' not in scrna_meta.columns:
        raise ValueError("❌ scRNA CSV must contain 'barcode' and 'celltype' columns")

    scrna_meta['barcode'] = scrna_meta['barcode'].astype(str).str.strip()
    scrna_meta['barcode_clean'] = scrna_meta['barcode'].str.replace('-1$', '', regex=True)
    scrna_meta.set_index('barcode_clean', inplace=True)

    # --- Prepare ATAC barcodes ---
    atac_barcodes = pd.Series(cistopic_obj.cell_data.index.astype(str), name='barcode')
    atac_barcodes_clean = atac_barcodes.str.replace('-1$', '', regex=True)
    cistopic_obj.cell_data['barcode'] = atac_barcodes.values
    cistopic_obj.cell_data['barcode_clean'] = atac_barcodes_clean.values
    cistopic_obj.cell_data.set_index('barcode_clean', inplace=True, drop=False)

    # --- Debug info ---
    print("\n=== BARCODE DEBUG INFO ===")
    print(f"ATAC total cells: {len(atac_barcodes)}")
    print(f"scRNA total cells: {len(scrna_meta)}")
    print(f"ATAC example barcodes: {atac_barcodes.head(5).tolist()}")
    print(f"scRNA example barcodes: {scrna_meta.index[:5].tolist()}")
    overlap = len(set(cistopic_obj.cell_data.index).intersection(scrna_meta.index))
    print(f"Shared barcodes (after cleaning): {overlap}")

    # --- Join scRNA metadata ---
    cistopic_obj.cell_data = cistopic_obj.cell_data.join(scrna_meta[['celltype']], how='left')

    # --- Report match statistics ---
    n_matched = cistopic_obj.cell_data['celltype'].notna().sum()
    total_cells = cistopic_obj.cell_data.shape[0]
    print(f"\n✅ Matched {n_matched} / {total_cells} cells ({n_matched / total_cells * 100:.2f}%)")

    unmatched = cistopic_obj.cell_data[cistopic_obj.cell_data['celltype'].isna()]
    if len(unmatched) > 0:
        print(f"⚠️ {len(unmatched)} ATAC cells without scRNA match. Showing 10 examples:")
        print(unmatched[['barcode', 'sample_id']].head(10).to_string(index=False))

    # --- Drop cells without metadata ---
    before_drop = cistopic_obj.cell_data.shape[0]
    cistopic_obj.cell_data.dropna(subset=['celltype'], inplace=True)
    after_drop = cistopic_obj.cell_data.shape[0]
    print(f"✓ Removed {before_drop - after_drop} cells without celltype info")

    # --- Remove empty celltype_atac column ---
    if 'celltype_atac' in cistopic_obj.cell_data.columns and cistopic_obj.cell_data['celltype_atac'].isna().all():
        cistopic_obj.cell_data.drop(columns=['celltype_atac'], inplace=True)
        print("✓ Removed empty 'celltype_atac' column")

    # --- Fix barcode suffixes like '-TH1___TH1' (optional) ---
    cistopic_obj.cell_data['barcode'] = cistopic_obj.cell_data['barcode'].apply(
        lambda x: re.sub(r'-TH\d+___TH\d+', '', str(x))
    )

    # --- Synchronize .cell_names length with filtered cells ---
    cistopic_obj.cell_names = cistopic_obj.cell_data.index.tolist()

    # --- Save cleaned object ---
    with open(output_pickle, "wb") as f:
        pickle.dump(cistopic_obj, f)

    print("\n💾 Updated and cleaned CistopicObject saved to:", output_pickle)
    print(f"📊 Final shape: {cistopic_obj.cell_data.shape}")
    print(f"📋 Columns: {list(cistopic_obj.cell_data.columns)}")
    print("✅ Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add scRNA metadata to merged CistopicObject (with cleanup)")
    parser.add_argument("--cistopic_pickle", required=True, help="Merged cisTopic pickle file")
    parser.add_argument("--scrna_csv", required=True, help="scRNA CSV with 'barcode' and 'celltype'")
    parser.add_argument("--output_pickle", required=True, help="Output path for updated pickle")
    args = parser.parse_args()

    add_scrna_metadata(args.cistopic_pickle, args.scrna_csv, args.output_pickle)

