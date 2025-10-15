#!/usr/bin/env python3
import pickle
import argparse
import pandas as pd
import os
import numpy as np
from scipy import sparse

def main():
    parser = argparse.ArgumentParser(description="Merge cisTopic objects while preserving original barcodes")
    parser.add_argument("input", help="Pickle file containing list/dict of cisTopic objects")
    parser.add_argument("output", help="Output path for merged cisTopic object")
    args = parser.parse_args()

    # Load input
    with open(args.input, "rb") as f:
        obj = pickle.load(f)

    # Convert to list if dict
    if isinstance(obj, dict):
        cistopic_list = list(obj.values())
    elif isinstance(obj, list):
        cistopic_list = obj
    else:
        raise TypeError(f"Expected list or dict of cisTopic objects, got {type(obj)}")

    print(f"✓ Processing {len(cistopic_list)} cisTopic objects")

    # Merge cell_data with original barcodes preserved
    merged_cell_data = []
    all_cell_names = []
    all_region_names = None
    
    for i, co in enumerate(cistopic_list):
        df = co.cell_data.copy()
        sample_name = f"Sample{i+1}"
        df["sample_id"] = sample_name
        
        # Store original barcode from index
        original_barcodes = df.index.tolist()
        df["orig_barcode"] = original_barcodes
        
        merged_cell_data.append(df)
        all_cell_names.extend(original_barcodes)
        print(f"✓ Sample {i+1}: {df.shape[0]} cells | Barcode examples: {original_barcodes[:3]}")

    # Concatenate while preserving original indices
    merged_df = pd.concat(merged_cell_data, axis=0)
    
    # Use original barcodes as index (they should be unique across samples)
    merged_df.index = merged_df["orig_barcode"]
    
    print(f"✓ Merged cell_data shape: {merged_df.shape}")

    # Merge region_data (should be identical across samples)
    if all(hasattr(co, 'region_data') for co in cistopic_list):
        first_region_data = cistopic_list[0].region_data
        if all(co.region_data.equals(first_region_data) for co in cistopic_list[1:]):
            merged_region_data = first_region_data.copy()
            all_region_names = merged_region_data.index.tolist()
            print(f"✓ Region data merged: {merged_region_data.shape}")
        else:
            print("⚠️  Region data differs between samples, using first sample's region_data")
            merged_region_data = first_region_data.copy()
            all_region_names = merged_region_data.index.tolist()
    else:
        # Get region names from counts matrix if region_data doesn't exist
        if hasattr(cistopic_list[0], 'counts'):
            all_region_names = cistopic_list[0].counts.columns.tolist()
            merged_region_data = pd.DataFrame(index=all_region_names)
        else:
            raise ValueError("No region data or counts matrix found to determine region names")

    # Merge count matrices
    merged_counts = None
    merged_binary = None
    merged_nr_frag = None
    merged_nr_acc = None
    
    # Merge counts matrices
    if all(hasattr(co, 'counts') and getattr(co, 'counts') is not None for co in cistopic_list):
        print("✓ Merging counts matrices...")
        count_matrices = []
        for co in cistopic_list:
            if sparse.issparse(co.counts):
                count_matrices.append(co.counts)
            else:
                count_matrices.append(sparse.csr_matrix(co.counts))
        merged_counts = sparse.vstack(count_matrices, format='csr')
        print(f"  Counts matrix merged: {merged_counts.shape}")

    # Merge binary matrices
    if all(hasattr(co, 'binary_count') and getattr(co, 'binary_count') is not None for co in cistopic_list):
        print("✓ Merging binary matrices...")
        binary_matrices = []
        for co in cistopic_list:
            if sparse.issparse(co.binary_count):
                binary_matrices.append(co.binary_count)
            else:
                binary_matrices.append(sparse.csr_matrix(co.binary_count))
        merged_binary = sparse.vstack(binary_matrices, format='csr')
        print(f"  Binary matrix merged: {merged_binary.shape}")

    # Merge nr_frag
    if all(hasattr(co, 'nr_frag') and getattr(co, 'nr_frag') is not None for co in cistopic_list):
        print("✓ Merging nr_frag...")
        nr_frag_list = []
        for co in cistopic_list:
            nr_frag_list.append(co.nr_frag)
        merged_nr_frag = np.concatenate(nr_frag_list)
        print(f"  nr_frag merged: {merged_nr_frag.shape}")

    # Merge nr_acc
    if all(hasattr(co, 'nr_acc') and getattr(co, 'nr_acc') is not None for co in cistopic_list):
        print("✓ Merging nr_acc...")
        nr_acc_list = []
        for co in cistopic_list:
            nr_acc_list.append(co.nr_acc)
        merged_nr_acc = np.concatenate(nr_acc_list)
        print(f"  nr_acc merged: {merged_nr_acc.shape}")

    # Handle path_to_fragments - combine if they exist
    merged_path_to_fragments = None
    fragment_paths = []
    for co in cistopic_list:
        if hasattr(co, 'path_to_fragments') and co.path_to_fragments is not None:
            fragment_paths.append(co.path_to_fragments)
    
    if fragment_paths:
        merged_path_to_fragments = fragment_paths[0] if len(fragment_paths) == 1 else fragment_paths
        print(f"✓ Fragment paths: {merged_path_to_fragments}")

    # Create new cisTopic object with required arguments
    try:
        from pycisTopic.cistopic_class import CistopicObject
        
        merged_obj = CistopicObject(
            fragment_matrix=merged_counts,
            binary_matrix=merged_binary,
            cell_names=all_cell_names,
            region_names=all_region_names,
            cell_data=merged_df,
            region_data=merged_region_data,
            path_to_fragments=merged_path_to_fragments
        )
        
        # Copy additional attributes if they exist
        for attr in ['project']:
            if hasattr(cistopic_list[0], attr):
                setattr(merged_obj, attr, getattr(cistopic_list[0], attr))
                
    except ImportError:
        # Fallback: use the first object's class directly
        merged_obj = cistopic_list[0].__class__(
            fragment_matrix=merged_counts,
            binary_matrix=merged_binary,
            cell_names=all_cell_names,
            region_names=all_region_names,
            cell_data=merged_df,
            region_data=merged_region_data,
            path_to_fragments=merged_path_to_fragments
        )

    # Save merged object
    out_path = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(out_path) if os.path.dirname(out_path) else '.', exist_ok=True)
    
    with open(out_path, "wb") as f:
        pickle.dump(merged_obj, f)

    print(f"\n✅ Merge completed successfully!")
    print(f"✓ Output: {out_path}")
    print(f"✓ Total cells: {merged_obj.cell_data.shape[0]}")
    print(f"✓ Cell_data columns: {list(merged_obj.cell_data.columns)}")
    
    # Show barcode overlap info
    if len(cistopic_list) == 2:
        barcodes1 = set(cistopic_list[0].cell_data.index)
        barcodes2 = set(cistopic_list[1].cell_data.index)
        shared = barcodes1.intersection(barcodes2)
        print(f"✓ Barcode overlap between samples: {len(shared)} shared barcodes")
        if len(shared) > 0:
            print(f"  Examples: {list(shared)[:5]}")
        print(f"  Unique to Sample1: {len(barcodes1 - barcodes2)}")
        print(f"  Unique to Sample2: {len(barcodes2 - barcodes1)}")

if __name__ == "__main__":
    main()
