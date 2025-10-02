#!/usr/bin/env python3
import os
import pickle
import argparse
import numpy as np
import pandas as pd
from pycisTopic.cistopic_class import CistopicObject

def debug_cistopic_object(pickle_path):
    """
    Debug a CistopicObject pickle file and print comprehensive information
    """
    print(f"=== DEBUGGING CISTOPIC OBJECT: {pickle_path} ===")
    print(f"File exists: {os.path.exists(pickle_path)}")
    
    if os.path.exists(pickle_path):
        file_size = os.path.getsize(pickle_path) / 1024  # Size in KB
        print(f"File size: {file_size:.2f} KB")
    
    # Load the pickle file
    try:
        with open(pickle_path, 'rb') as f:
            data = pickle.load(f)
        print(f"\n✅ Successfully loaded pickle file")
    except Exception as e:
        print(f"❌ Failed to load pickle: {e}")
        return
    
    # Determine data type
    print(f"\n=== DATA TYPE: {type(data)} ===")
    
    if isinstance(data, list):
        print(f"List length: {len(data)}")
        for i, item in enumerate(data):
            print(f"  [{i}] Type: {type(item)}")
            if hasattr(item, '__class__'):
                print(f"      Class: {item.__class__.__name__}")
        
        # If it's a list of CistopicObjects, debug each one
        if all(isinstance(item, CistopicObject) for item in data):
            print(f"\n=== DETAILED CISTOPIC OBJECT ANALYSIS ===")
            for i, cistopic_obj in enumerate(data):
                print(f"\n--- CistopicObject [{i}] ---")
                debug_single_cistopic(cistopic_obj, f"Sample_{i}")
        
    elif isinstance(data, CistopicObject):
        print(f"\n=== SINGLE CISTOPIC OBJECT ANALYSIS ===")
        debug_single_cistopic(data, "Single_Object")
    else:
        print(f"Unexpected data type: {type(data)}")
        print(f"Available attributes: {dir(data)}")

def debug_single_cistopic(cistopic_obj, obj_name):
    """
    Debug a single CistopicObject instance
    """
    print(f"\n📊 Debugging: {obj_name}")
    
    # Basic attributes
    print(f"  Project name: {getattr(cistopic_obj, 'project', 'Not available')}")
    print(f"  Path to fragments: {getattr(cistopic_obj, 'path_to_fragments', 'Not available')}")
    
    # Check main data components
    components = ['cell_data', 'region_data', 'count_matrix', 'binary_count_matrix', 
                  'projection_data', 'topic_binarized_matrix', 'topic_region']
    
    for component in components:
        if hasattr(cistopic_obj, component):
            attr_value = getattr(cistopic_obj, component)
            if attr_value is not None:
                print(f"\n  🔍 {component}:")
                if hasattr(attr_value, 'shape'):
                    print(f"      Shape: {attr_value.shape}")
                if hasattr(attr_value, 'columns'):
                    print(f"      Columns: {len(attr_value.columns)}")
                    if len(attr_value.columns) > 0:
                        print(f"      First 5 columns: {list(attr_value.columns[:5])}")
                if hasattr(attr_value, 'index'):
                    print(f"      Index length: {len(attr_value.index)}")
                    if len(attr_value.index) > 0:
                        print(f"      First 5 index: {list(attr_value.index[:5])}")
                
                # For sparse matrices
                if hasattr(attr_value, 'nnz'):
                    print(f"      Non-zero elements: {attr_value.nnz}")
                    if hasattr(attr_value, 'shape'):
                        density = attr_value.nnz / (attr_value.shape[0] * attr_value.shape[1])
                        print(f"      Density: {density:.6f}")
                
                # For dataframes, show some statistics
                if isinstance(attr_value, pd.DataFrame):
                    print(f"      Data types: {attr_value.dtypes.value_counts().to_dict()}")
                    if not attr_value.empty:
                        print(f"      First few rows:")
                        print(attr_value.head(3))
            else:
                print(f"  🔍 {component}: None")
        else:
            print(f"  🔍 {component}: Not present in object")
    
    # Check for topic model data
    topic_attrs = ['cell_topic', 'topic_region', 'topic_binarized_matrix']
    for attr in topic_attrs:
        if hasattr(cistopic_obj, attr) and getattr(cistopic_obj, attr) is not None:
            print(f"\n  🎯 Topic Model - {attr}:")
            attr_value = getattr(cistopic_obj, attr)
            if hasattr(attr_value, 'shape'):
                print(f"      Shape: {attr_value.shape}")
            if isinstance(attr_value, dict):
                print(f"      Dictionary keys: {list(attr_value.keys())}")
                for key, value in attr_value.items():
                    if hasattr(value, 'shape'):
                        print(f"        {key} shape: {value.shape}")
    
    # Check for dimensionality reduction
    if hasattr(cistopic_obj, 'projection_data') and cistopic_obj.projection_data is not None:
        print(f"\n  📈 Projection Data:")
        for method, projection in cistopic_obj.projection_data.items():
            print(f"      {method}: {projection.shape if hasattr(projection, 'shape') else 'No shape'}")
    
    # Summary statistics
    print(f"\n  📋 SUMMARY STATISTICS:")
    if hasattr(cistopic_obj, 'cell_data') and cistopic_obj.cell_data is not None:
        print(f"      Number of cells: {len(cistopic_obj.cell_data)}")
    if hasattr(cistopic_obj, 'region_data') and cistopic_obj.region_data is not None:
        print(f"      Number of regions: {len(cistopic_obj.region_data)}")
    if hasattr(cistopic_obj, 'count_matrix') and cistopic_obj.count_matrix is not None:
        print(f"      Count matrix shape: {cistopic_obj.count_matrix.shape}")
        if hasattr(cistopic_obj.count_matrix, 'nnz'):
            print(f"      Total counts: {cistopic_obj.count_matrix.sum()}")

def main():
    parser = argparse.ArgumentParser(description='Debug CistopicObject pickle files')
    parser.add_argument('input', help='Path to CistopicObject pickle file')
    args = parser.parse_args()
    
    debug_cistopic_object(args.input)

if __name__ == "__main__":
    main()
