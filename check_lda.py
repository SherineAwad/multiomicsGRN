#!/usr/bin/env python3
import argparse
import pickle
import numpy as np
import pandas as pd

def check_lda_model(pickle_path):
    print("✅ Loading CistopicObject...")
    with open(pickle_path, "rb") as f:
        obj = pickle.load(f)

    print("\n=== BASIC INFO ===")
    print("Data type:", type(obj))
    print("Project:", getattr(obj, "project", "N/A"))

    print("\n=== MODEL PRESENCE ===")
    models = getattr(obj, "models", None)
    selected_model = getattr(obj, "selected_model", None)
    print("Raw type of obj.models:", type(models))
    print("Selected model:", selected_model)

    if models is None:
        print("❌ No 'models' attribute found.")
        return

    if isinstance(models, dict):
        print("Available model keys:", list(models.keys()) if models else "❌ Empty dict")
    elif isinstance(models, list):
        print(f"⚠️ 'models' is a list with length {len(models)}")
        for i, m in enumerate(models):
            print(f"  [{i}] Type: {type(m)}")
    else:
        print(f"⚠️ 'models' is of unexpected type: {type(models)}")

    # topic_region
    tr = getattr(obj, "topic_region", None)
    if tr is None:
        print("\n❌ No topic_region found. LDA not properly added.")
    else:
        print("\n=== TOPIC_REGION MATRIX ===")
        print("Shape:", tr.shape)
        print("Max:", tr.max().max(), "Min:", tr.min().min())
        print("Total nonzero entries:", (tr.values > 0).sum(), "of", tr.size)
        print("Values > 1e-4:", (tr.values > 1e-4).sum())

    # cell_topic
    ct = getattr(obj, "cell_topic", None)
    if ct is None:
        print("\n❌ No cell_topic found.")
    else:
        print("\n=== CELL_TOPIC MATRIX ===")
        print("Shape:", ct.shape)
        print("NaNs:", np.isnan(ct.values).sum())
        print("Per-topic mean (first 10):")
        print(ct.mean(axis=1).round(6).head(10))

    print("\n✅ Check complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check if LDA model was added correctly to CistopicObject before binarisation.")
    parser.add_argument("pickle_path", help="Path to CistopicObject pickle file")
    args = parser.parse_args()
    check_lda_model(args.pickle_path)

