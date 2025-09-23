#!/usr/bin/env python3
import os
import pickle
import argparse

def main():
    parser = argparse.ArgumentParser(description="Add the existing selected_model to the same CistopicObject")
    parser.add_argument("--input_pickle", "-i", required=True, help="Path to CistopicObject pickle file")
    parser.add_argument("--output_dir", "-o", required=True, help="Directory to save updated object")
    args = parser.parse_args()

    # Load the CistopicObject
    with open(args.input_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    # Add the existing selected_model to the object
    if hasattr(cistopic_obj, "selected_model") and cistopic_obj.selected_model is not None:
        cistopic_obj.add_LDA_model(cistopic_obj.selected_model)
        print("Selected LDA model added to object.")
    else:
        print("[WARNING] No selected_model found in object. Skipping add_LDA_model step.")

    # Check if project method exists
    if hasattr(cistopic_obj, "project") and callable(getattr(cistopic_obj, "project", None)):
        # Safe to call project
        cistopic_obj.project(target="cell")
        print("Projection to cells done successfully.")
    else:
        print("[WARNING] The 'project' method is missing or overwritten. Cannot run cistopic_obj.project().")

    # Save updated object
    os.makedirs(args.output_dir, exist_ok=True)
    out_path = os.path.join(args.output_dir, "cistopic_LDA_obj.pkl")
    with open(out_path, "wb") as f:
        pickle.dump(cistopic_obj, f)

    print(f"Updated CistopicObject saved to: {out_path}")

if __name__ == "__main__":
    main()

