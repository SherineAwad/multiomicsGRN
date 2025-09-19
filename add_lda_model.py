#!/usr/bin/env python3
import os
import pickle
import argparse

def main():
    parser = argparse.ArgumentParser(description="Add the existing selected_model to the same CistopicObject")
    parser.add_argument("--cistopic_obj", "-c", required=True, help="Path to CistopicObject pickle file")
    parser.add_argument("--out_dir", "-o", required=True, help="Directory to save updated object")
    args = parser.parse_args()

    # Load the CistopicObject
    with open(args.cistopic_obj, "rb") as f:
        cistopic_obj = pickle.load(f)

    # Add the existing selected_model to the object
    cistopic_obj.add_LDA_model(cistopic_obj.selected_model)

    # Save updated object
    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, "cistopic_obj.pkl")
    with open(out_path, "wb") as f:
        pickle.dump(cistopic_obj, f)

    print(f"Existing selected_model added again. Updated object saved to: {out_path}")

if __name__ == "__main__":
    main()

