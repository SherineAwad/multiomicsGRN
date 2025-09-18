#!/usr/bin/env python3
import argparse
import pickle
import os
from pycisTopic.cistopic_class import merge

def main():
    parser = argparse.ArgumentParser(
        description="Merge cisTopic objects stored inside a single pickle file."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to pickle file containing multiple cisTopic objects."
    )
    parser.add_argument(
        "--output",
        required=True,
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

    # Merge
    print("Merging cisTopic objects...")
    merged_obj = merge(cistopic_list)

    # Save merged
    out_path = os.path.abspath(args.output)
    with open(out_path, "wb") as fh:
        pickle.dump(merged_obj, fh)
    print(f"Merged cisTopic object saved to: {out_path}")

if __name__ == "__main__":
    main()

