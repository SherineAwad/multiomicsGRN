#!/usr/bin/env python3
import pickle
import argparse

def inspect_qc_pickle(pickle_path):
    """
    Inspect SCENIC+ QC barcodes pickle file
    """
    print(f"=== Inspecting QC pickle: {pickle_path} ===")
    
    with open(pickle_path, 'rb') as f:
        data = pickle.load(f)
    
    print(f"Top-level keys: {list(data.keys())}")
    print()
    
    print("=== Number of barcodes per sample ===")
    barcodes_dict = data['barcodes']
    for sample, barcodes in barcodes_dict.items():
        print(f"{sample}: {len(barcodes)} barcodes")
    
    print()
    print("=== QC thresholds per sample ===")
    thresholds_dict = data['thresholds']
    for sample, thresholds in thresholds_dict.items():
        print(f"{sample}: {thresholds}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Inspect QC barcodes pickle")
    parser.add_argument("pickle_file", help="Path to QC pickle file")
    args = parser.parse_args()
    
    inspect_qc_pickle(args.pickle_file)
