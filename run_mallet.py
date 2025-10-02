#!/usr/bin/env python3
import os
import argparse
import pickle
from pycisTopic.lda_models import run_cgs_models_mallet
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments

def run_mallet_on_object(cistopic_obj, mallet_path, n_topics, n_cpu, n_iter,
                         tmp_path, save_path, mallet_memory, random_state,
                         alpha, alpha_by_topic, eta, eta_by_topic):
    """
    Run MALLET on a single merged CistopicObject and select the best model.
    """
    sample_name = "merged_sample"
    print(f"Running MALLET for {sample_name}...")

    # Ensure tmp and save directories exist
    sample_tmp_path = os.path.join(tmp_path, sample_name)
    sample_save_path = os.path.join(save_path, sample_name)
    os.makedirs(sample_tmp_path, exist_ok=True)
    os.makedirs(sample_save_path, exist_ok=True)

    # Set Java memory for Mallet
    os.environ['MALLET_MEMORY'] = mallet_memory

    # Run MALLET CGS models
    models = run_cgs_models_mallet(
        cistopic_obj,
        n_topics=n_topics,
        n_cpu=n_cpu,
        n_iter=n_iter,
        random_state=random_state,
        alpha=alpha,
        alpha_by_topic=alpha_by_topic,
        eta=eta,
        eta_by_topic=eta_by_topic,
        tmp_path=sample_tmp_path,
        save_path=sample_save_path,
        mallet_path=mallet_path,
    )

    if not models or len(models) == 0:
        raise ValueError("No MALLET models were generated.")

    # Select best model by log likelihood
    best_model = max(models, key=lambda m: getattr(m, 'log_likelihood', -1e10))

    # Add model to CistopicObject
    cistopic_obj.add_LDA_model(best_model)
    cistopic_obj.selected_model = best_model

    print(f"Finished MALLET. Model added to merged CistopicObject.\n")
    return cistopic_obj

def main(cistopic_obj_pickle, mallet_path, n_topics, n_cpu, n_iter,
         tmp_path, save_path, mallet_memory="250G", random_state=555,
         alpha=50, alpha_by_topic=True, eta=0.1, eta_by_topic=False):

    # Load merged CistopicObject
    with open(cistopic_obj_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    # FIX: Check if count_matrix exists using hasattr instead
    if not hasattr(cistopic_obj, 'count_matrix') or cistopic_obj.count_matrix is None:
        print("Creating count matrix for Mallet...")
        
        # Extract regions from region_data
        regions_list = cistopic_obj.region_data.index.tolist()
        
        # Get cell barcodes
        cell_barcodes = cistopic_obj.cell_data.index.tolist()
        
        # FIX: Use the actual fragment file paths from the dictionary
        # Since we have multiple samples, we need to handle this differently
        # Let's use the first fragment file for simplicity, or create a combined approach
        fragment_files = list(cistopic_obj.path_to_fragments.values())
        if len(fragment_files) > 1:
            print(f"Warning: Multiple fragment files found. Using first one: {fragment_files[0]}")
        
        # Create a new CistopicObject with count matrix using the first fragment file
        cistopic_obj_with_matrix = create_cistopic_object_from_fragments(
            path_to_fragments=fragment_files[0],  # Use string path, not dict
            path_to_regions=regions_list,
            metrics=cistopic_obj.cell_data,
            valid_bc=cell_barcodes,
            project=cistopic_obj.project,
            split_pattern='-',
            n_cpu=n_cpu
        )
        
        # Replace the original object with the one that has count matrix
        cistopic_obj = cistopic_obj_with_matrix
        print("✓ Count matrix created successfully")
    else:
        print("✓ Count matrix already exists")

    # Run MALLET
    updated_obj = run_mallet_on_object(
        cistopic_obj=cistopic_obj,
        mallet_path=mallet_path,
        n_topics=n_topics,
        n_cpu=n_cpu,
        n_iter=n_iter,
        tmp_path=tmp_path,
        save_path=save_path,
        mallet_memory=mallet_memory,
        random_state=random_state,
        alpha=alpha,
        alpha_by_topic=alpha_by_topic,
        eta=eta,
        eta_by_topic=eta_by_topic
    )

    # Save updated object
    output_pickle = os.path.join(save_path, "merged_cistopic_with_models.pkl")
    with open(output_pickle, "wb") as f:
        pickle.dump(updated_obj, f)
    print(f"✓ Updated CistopicObject with models saved to: {output_pickle}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MALLET topic modeling on a merged CistopicObject")
    parser.add_argument("--cistopic_obj_pickle", required=True)
    parser.add_argument("--mallet_path", required=True)
    parser.add_argument("--n_topics", type=int, nargs="+", required=True)
    parser.add_argument("--n_cpu", type=int, default=12)
    parser.add_argument("--n_iter", type=int, default=500)
    parser.add_argument("--tmp_path", required=True)
    parser.add_argument("--save_path", required=True)
    parser.add_argument("--mallet_memory", default="300G")
    parser.add_argument("--random_state", type=int, default=555)
    parser.add_argument("--alpha", type=float, default=5.0)
    parser.add_argument("--alpha_by_topic", action="store_true", default=True)
    parser.add_argument("--eta", type=float, default=0.1)
    parser.add_argument("--eta_by_topic", action="store_true", default=True)
    args = parser.parse_args()

    main(
        cistopic_obj_pickle=args.cistopic_obj_pickle,
        mallet_path=args.mallet_path,
        n_topics=args.n_topics,
        n_cpu=args.n_cpu,
        n_iter=args.n_iter,
        tmp_path=args.tmp_path,
        save_path=args.save_path,
        mallet_memory=args.mallet_memory,
        random_state=args.random_state,
        alpha=args.alpha,
        alpha_by_topic=args.alpha_by_topic,
        eta=args.eta,
        eta_by_topic=args.eta_by_topic
    )
