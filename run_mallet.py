#!/usr/bin/env python3
import os
import argparse
import pickle
from pycisTopic.lda_models import run_cgs_models_mallet

def run_mallet_on_object(cistopic_obj, mallet_path, n_topics, n_cpu, n_iter,
                         tmp_path, save_path, mallet_memory, random_state,
                         alpha, alpha_by_topic, eta, eta_by_topic, sample_name):
    """
    Run MALLET on a single CistopicObject and select the best model.
    """
    print(f"Running MALLET for sample: {sample_name}")

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

    # Select the best model automatically
    if not models:
        raise ValueError(f"No models were generated for sample {sample_name}")
    best_model = models[0]  # pycisTopic usually selects by highest likelihood internally
    cistopic_obj.add_LDA_model(best_model)
    cistopic_obj.selected_model = best_model

    print(f"Finished MALLET for {sample_name}. Model added to CistopicObject.\n")
    return cistopic_obj

def main(cistopic_obj_pickle, mallet_path, n_topics, n_cpu, n_iter,
         tmp_path, save_path, mallet_memory="250G", random_state=555,
         alpha=50, alpha_by_topic=True, eta=0.1, eta_by_topic=False):
    """
    Run Mallet-based topic modeling on a CistopicObject or list of objects.
    """
    # Load CistopicObject(s)
    if not os.path.exists(cistopic_obj_pickle):
        raise FileNotFoundError(f"CistopicObject pickle not found: {cistopic_obj_pickle}")
    with open(cistopic_obj_pickle, "rb") as f:
        cistopic_objs = pickle.load(f)

    # Ensure we have a list
    if not isinstance(cistopic_objs, list):
        cistopic_objs = [cistopic_objs]

    updated_objs = []
    for idx, obj in enumerate(cistopic_objs):
        sample_name = f"sample_{idx}"
        updated_obj = run_mallet_on_object(
            cistopic_obj=obj,
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
            eta_by_topic=eta_by_topic,
            sample_name=sample_name
        )
        updated_objs.append(updated_obj)

    # Save updated list back to pickle
    output_pickle = os.path.join(save_path, "cistopic_objs_with_models.pkl")
    with open(output_pickle, "wb") as f:
        pickle.dump(updated_objs, f)

    print(f"All samples processed. Updated pickle saved to: {output_pickle}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MALLET topic modeling on a CistopicObject or list of objects")
    parser.add_argument("--cistopic_obj_pickle", required=True,
                        help="Pickle file containing a CistopicObject or list of objects")
    parser.add_argument("--mallet_path", required=True,
                        help="Full path to Mallet executable (bin/mallet)")
    parser.add_argument("--n_topics", type=int, nargs="+",
                        default=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
                        help="List of topic numbers to run")
    parser.add_argument("--n_cpu", type=int, default=12, help="Number of CPUs to use")
    parser.add_argument("--n_iter", type=int, default=500, help="Number of iterations for Mallet")
    parser.add_argument("--tmp_path", required=True, help="Temporary folder for intermediate files")
    parser.add_argument("--save_path", required=True, help="Folder to save resulting models")
    parser.add_argument("--mallet_memory", default="200G", help="Max Java memory for Mallet")
    parser.add_argument("--random_state", type=int, default=555, help="Random seed")
    parser.add_argument("--alpha", type=float, default=50, help="Alpha hyperparameter")
    parser.add_argument("--alpha_by_topic", type=bool, default=True, help="Use alpha per topic")
    parser.add_argument("--eta", type=float, default=0.1, help="Eta hyperparameter")
    parser.add_argument("--eta_by_topic", type=bool, default=False, help="Use eta per topic")
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

