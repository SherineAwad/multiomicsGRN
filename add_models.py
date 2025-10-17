#!/usr/bin/env python3
import os
import pickle
import glob
import argparse
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='Add LDA models to CistopicObject')
    parser.add_argument('-i', '--input', required=True, help='Input CistopicObject pickle')
    parser.add_argument('-o', '--output', required=True, help='Output CistopicObject pickle')
    args = parser.parse_args()

    # Load CistopicObject
    with open(args.input, 'rb') as f:
        cistopic_obj = pickle.load(f)

    # Find all topic model files
    mallet_dir = os.path.dirname(args.input)
    model_files = glob.glob(os.path.join(mallet_dir, "**/Topic*.pkl"), recursive=True)

    if not model_files:
        raise FileNotFoundError("No Topic*.pkl files found")

    print(f"Found {len(model_files)} topic models")

    # Load all models
    models = []
    for model_file in model_files:
        with open(model_file, 'rb') as f:
            model = pickle.load(f)
            models.append(model)
            print(f"Loaded: {os.path.basename(model_file)}")

    # Select best model
    best_model = max(models, key=lambda m: m.metrics.loc['Metric', 'loglikelihood'])

    # Manually attach the topic matrices - TRANSPOSE cell_topic to (cells, topics)
    cistopic_obj.cell_topic = best_model.cell_topic.T
    cistopic_obj.topic_region = best_model.topic_region
    cistopic_obj.models = models
    cistopic_obj.selected_model = best_model

    print(f"✓ Manually attached topic matrices from {best_model.n_topic} topic model")
    print(f"  cell_topic shape: {cistopic_obj.cell_topic.shape}")
    print(f"  topic_region shape: {cistopic_obj.topic_region.shape}")

    # Save updated object
    with open(args.output, 'wb') as f:
        pickle.dump(cistopic_obj, f)
    print(f"Saved to: {args.output}")

if __name__ == "__main__":
    main()
