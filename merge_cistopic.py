#!/usr/bin/env python3
import pickle
import argparse
import os
import time
from scipy import sparse
import pandas as pd
import numpy as np
from pycisTopic.cistopic_class import CistopicObject

def main():
    parser = argparse.ArgumentParser(description="Merge cisTopic objects while preserving original barcodes")
    parser.add_argument("input", help="Pickle file containing list/dict of cisTopic objects")
    parser.add_argument("output", help="Output path for merged cisTopic object")
    args = parser.parse_args()

    # Load input
    print("Loading cisTopic objects...")
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

    # Quick dimension check
    for i, co in enumerate(cistopic_list):
        print(f"  Object {i+1}: {len(co.cell_names)} cells × {len(co.region_names)} regions")

    # Merge cell_data
    print("Merging cell_data...")
    merged_cell_data = []
    all_cell_names = []
    for i, co in enumerate(cistopic_list):
        df = co.cell_data.copy()
        df["sample_id"] = f"Sample{i+1}"
        df["orig_barcode"] = df.index.tolist()
        merged_cell_data.append(df)
        all_cell_names.extend(df.index.tolist())
    merged_cell_data = pd.concat(merged_cell_data, axis=0)
    print(f"✓ Merged cell_data: {merged_cell_data.shape}")

    # Collect all unique regions
    print("Collecting unique regions...")
    all_regions = sorted(set().union(*[co.region_data.index.tolist() for co in cistopic_list]))
    region_index_map = {r: i for i, r in enumerate(all_regions)}
    print(f"✓ Total unique regions: {len(all_regions)}")

    # FAST VERSION: Merge matrices using efficient alignment
    print("Merging matrices (fast method)...")
    
    def merge_matrix_fast(matrix_name):
        """Fast matrix merging using column alignment"""
        aligned_matrices = []
        
        for obj_idx, co in enumerate(cistopic_list):
            print(f"  Processing {matrix_name} for object {obj_idx+1}...")
            start_time = time.time()
            
            # Get the matrix
            if matrix_name == "fragment_matrix":
                mat = getattr(co, "fragment_matrix", None)
            else:  # binary_matrix
                mat = getattr(co, "binary_matrix", None)
                if mat is None:
                    # Create binary from fragment_matrix if needed
                    mat_source = getattr(co, "fragment_matrix", None)
                    if mat_source is not None:
                        mat = (mat_source > 0).astype(int)
            
            if mat is None:
                raise ValueError(f"CistopicObject has no {matrix_name}: {co}")
            
            mat = sparse.csr_matrix(mat)
            print(f"    Original shape: {mat.shape}")
            
            # Create mapping from object regions to merged regions
            obj_region_to_idx = {region: idx for idx, region in enumerate(co.region_data.index)}
            
            # Find which regions from this object exist in the merged set
            valid_regions = [r for r in co.region_data.index if r in region_index_map]
            print(f"    Found {len(valid_regions)}/{len(co.region_data.index)} regions in merged set")
            
            if not valid_regions:
                print("    WARNING: No regions match between objects!")
                # Create empty matrix as fallback
                new_mat = sparse.csr_matrix((mat.shape[0], len(all_regions)))
            else:
                # Get column indices in both original and merged matrices
                obj_col_indices = [obj_region_to_idx[r] for r in valid_regions]
                merged_col_indices = [region_index_map[r] for r in valid_regions]
                
                # Filter to valid column indices
                valid_obj_indices = [idx for idx in obj_col_indices if idx < mat.shape[1]]
                valid_merged_indices = [merged_col_indices[i] for i, idx in enumerate(obj_col_indices) if idx < mat.shape[1]]
                
                # Create new matrix efficiently
                new_mat = sparse.lil_matrix((mat.shape[0], len(all_regions)))
                
                # Copy valid columns in batch
                if valid_obj_indices:
                    new_mat[:, valid_merged_indices] = mat[:, valid_obj_indices]
            
            aligned_matrices.append(new_mat.tocsr())
            print(f"    ✓ Finished in {time.time() - start_time:.1f}s")
        
        # Stack all matrices
        merged_matrix = sparse.vstack(aligned_matrices, format='csr')
        print(f"✓ Merged {matrix_name} shape: {merged_matrix.shape}")
        return merged_matrix

    # Merge both matrices using fast method
    merged_counts = merge_matrix_fast("fragment_matrix")
    merged_binary = merge_matrix_fast("binary_matrix")

    # Merge region_data
    print("Merging region_data...")
    merged_region_data = pd.DataFrame(index=all_regions)

    # Merge path_to_fragments
    merged_path_to_fragments = []
    for co in cistopic_list:
        if getattr(co, "path_to_fragments", None) is not None:
            merged_path_to_fragments.append(co.path_to_fragments)
    if len(merged_path_to_fragments) == 1:
        merged_path_to_fragments = merged_path_to_fragments[0]

    # Create new CistopicObject with matrices
    print("Creating merged cisTopic object...")
    merged_obj = CistopicObject(
        fragment_matrix=merged_counts,
        binary_matrix=merged_binary,
        cell_names=all_cell_names,
        region_names=all_regions,
        cell_data=merged_cell_data,
        region_data=merged_region_data,
        path_to_fragments=merged_path_to_fragments
    )

    # Save merged object
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "wb") as f:
        pickle.dump(merged_obj, f)

    print(f"\n✅ Merge completed. Output saved to {args.output}")
    print(f"✓ Total cells: {merged_obj.cell_data.shape[0]}")
    print(f"✓ Total regions: {len(merged_obj.region_names)}")
    print(f"✓ fragment_matrix shape: {merged_obj.fragment_matrix.shape}")
    print(f"✓ binary_matrix shape: {merged_obj.binary_matrix.shape}")

if __name__ == "__main__":
    main()
