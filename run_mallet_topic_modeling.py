#!/usr/bin/env python3
import os
import argparse
import pickle
from pycisTopic.lda_models import run_cgs_models_mallet

def main(
    cistopic_obj_pickle,
    mallet_path,
    n_topics,
    n_cpu,
    n_iter,
    tmp_path,
    save_path,
    mallet_memory="200G",
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False
):
    """
    Run Mallet-based topic modeling on a CistopicObject.
    
    cistopic_obj_pickle: path to CistopicObject pickle
    mallet_path: full path to Mallet executable
    n_topics: list of topic numbers to test (e.g., [2,5,10])
    n_cpu: number of CPUs to use
    n_iter: number of iterations
    tmp_path: temporary folder for intermediate files
    save_path: folder to save resulting models
    mallet_memory: max memory for Mallet (default '200G')
    """

    # Set Java memory for Mallet
    os.environ['MALLET_MEMORY'] = mallet_memory

    # Load CistopicObject
    if not os.path.exists(cistopic_obj_pickle):
        raise FileNotFoundError(f"CistopicObject pickle not found: {cistopic_obj_pickle}")
    with open(cistopic_obj_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    # Ensure tmp and save directories exist
    os.makedirs(tmp_path, exist_ok=True)
    os.makedirs(save_path, exist_ok=True)

    # Run Mallet CGS models
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
        tmp_path=tmp_path,
        save_path=save_path,
        mallet_path=mallet_path,
    )

    print(f"Finished running Mallet topic models. Models saved in: {save_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Mallet-based CGS topic modeling on a CistopicObject")
    parser.add_argument("--cistopic_obj_pickle", required=True,
                        help="Pickle file containing a CistopicObject")
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

