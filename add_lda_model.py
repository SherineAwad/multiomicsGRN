#!/usr/bin/env python3
import os
import pickle
import argparse
from pycisTopic.lda_models import evaluate_models  # import evaluate_models

def main():
    parser = argparse.ArgumentParser(description="Add the existing selected_model to the same CistopicObject")
    parser.add_argument("--input_pickle", "-i", required=True, help="Path to CistopicObject pickle file")
    parser.add_argument("--output_dir", "-o", required=True, help="Directory to save updated object")
    parser.add_argument("--select_model", "-s", type=int, default=None,
                        help="Number of topics to select after evaluation (optional)")
    args = parser.parse_args()

    # Load the CistopicObject
    with open(args.input_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    # Evaluate models if available
    if hasattr(cistopic_obj, "models") and cistopic_obj.models:
        print("Evaluating LDA models...")
        selected_model = evaluate_models(
            cistopic_obj.models,
            select_model=args.select_model,
            return_model=True
        )
        # Save selected model to the object
        cistopic_obj.selected_model = selected_model
        print(f"Selected LDA model (topics={getattr(selected_model, 'n_topics', getattr(selected_model, 'n_topic', 'N/A'))}) stored in object.")
    else:
        print("[WARNING] No models found in CistopicObject to evaluate.")
        # Fallback: keep existing selected_model if present
        if not hasattr(cistopic_obj, "selected_model") or cistopic_obj.selected_model is None:
            print("[WARNING] No selected_model found. Please check your object.")

    # Add the selected_model to the object safely
    if hasattr(cistopic_obj, "selected_model") and cistopic_obj.selected_model is not None:
        try:
            cistopic_obj.add_LDA_model(cistopic_obj.selected_model)
            print("Selected LDA model added to object.")
        except Exception as e:
            print(f"[WARNING] Could not add LDA model: {e}")

    # Attempt projection if method exists
    project_method = getattr(cistopic_obj, "project", None)
    if callable(project_method):
        try:
            cistopic_obj.project(target="cell")
            print("Projection to cells done successfully.")
        except Exception as e:
            print(f"[WARNING] 'project' method failed: {e}")
    else:
        print("[INFO] 'project' method is missing or not callable. Skipping projection step.")

    # Save updated object
    os.makedirs(args.output_dir, exist_ok=True)
    out_path = os.path.join(args.output_dir, "cistopic_LDA_obj.pkl")
    with open(out_path, "wb") as f:
        pickle.dump(cistopic_obj, f)

    print(f"Updated CistopicObject saved to: {out_path}")

if __name__ == "__main__":
    main()

