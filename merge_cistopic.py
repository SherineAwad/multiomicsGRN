#!/usr/bin/env python3
import argparse
import pickle
import os
import copy

def main():
    parser = argparse.ArgumentParser(
        description="Merge cisTopic objects stored inside a single pickle file."
    )
    parser.add_argument(
        "input",
        help="Path to pickle file containing multiple cisTopic objects."
    )
    parser.add_argument(
        "output", 
        help="Path to save merged cisTopic object (e.g., merged_cistopic.pkl)."
    )
    args = parser.parse_args()

    # Load pickle
    print(f"Loading cisTopic objects from {args.input}...")
    with open(args.input, "rb") as fh:
        obj = pickle.load(fh)

    # obj might be a list or dict — handle both
    if isinstance(obj, dict):
        cistopic_list = list(obj.values())
        print(f"Found {len(cistopic_list)} cisTopic objects (dict keys: {list(obj.keys())})")
    elif isinstance(obj, list):
        cistopic_list = obj
        print(f"Found {len(cistopic_list)} cisTopic objects (list)")
    else:
        raise TypeError(
            f"Unexpected object type in pickle: {type(obj)}. Expected list or dict of cisTopic objects."
        )

    # Method 1: Try with copy=True first
    print("Merging cisTopic objects...")
    try:
        merged_obj = cistopic_list[0].merge(
            cistopic_obj_list=cistopic_list[1:],
            project='TH1_TH2_merged',
            split_pattern='-',
            copy=True  # This ensures the method returns a new object
        )
        print("✓ Merge completed with copy=True")
    except Exception as e:
        print(f"Merge with copy=True failed: {e}")
        print("Trying alternative method with deep copy...")
        
        # Method 2: Create a deep copy first and merge in-place
        base_obj = copy.deepcopy(cistopic_list[0])
        merged_obj = base_obj.merge(
            cistopic_obj_list=cistopic_list[1:],
            project='TH1_TH2_merged', 
            split_pattern='-',
            copy=False  # Modify the copy in-place
        )
        print("✓ Merge completed with deep copy method")

    # Verify the merged object is not None before saving
    if merged_obj is None:
        raise ValueError("Merge operation returned None. Check the merge method parameters.")

    # Save merged object
    out_path = os.path.abspath(args.output)
    with open(out_path, "wb") as fh:
        pickle.dump(merged_obj, fh)
    
    # Verify the saved file
    if os.path.exists(out_path):
        file_size = os.path.getsize(out_path) / 1024  # Size in KB
        print(f"✓ Merged cisTopic object saved to: {out_path}")
        print(f"✓ File size: {file_size:.2f} KB")
        
        # Quick verification load
        with open(out_path, "rb") as fh:
            verify_obj = pickle.load(fh)
        print(f"✓ Verification: Loaded object type: {type(verify_obj)}")
    else:
        print("✗ Error: Output file was not created")

if __name__ == "__main__":
    main()
