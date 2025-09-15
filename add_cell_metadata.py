#!/usr/bin/env python3
import os
import argparse
import pickle
import pandas as pd

def main(cistopic_objects_pickle, cell_data_tsv, output_dir, split_pattern="-"):
    """
    Add cell metadata to CistopicObjects and save them back to pickle files.
    
    cistopic_objects_pickle: pickle file containing a list of CistopicObjects
    cell_data_tsv: TSV file with cell metadata (rows = barcodes, columns = metadata)
    output_dir: directory to save annotated CistopicObjects
    split_pattern: pattern to split barcodes (default '-')
    """

    # Load cistopic objects
    if not os.path.exists(cistopic_objects_pickle):
        raise FileNotFoundError(f"cistopic_objects pickle not found: {cistopic_objects_pickle}")
    with open(cistopic_objects_pickle, "rb") as f:
        cistopic_obj_list = pickle.load(f)

    # Load cell metadata
    if not os.path.exists(cell_data_tsv):
        raise FileNotFoundError(f"Cell metadata file not found: {cell_data_tsv}")
    cell_data = pd.read_table(cell_data_tsv, index_col=0)

    # Create output dir
    os.makedirs(output_dir, exist_ok=True)

    # Add cell data to each CistopicObject
    for cistopic_obj in cistopic_obj_list:
        sample_id = cistopic_obj.project
        print(f"Adding metadata to sample: {sample_id}")

        # Add metadata
        cistopic_obj.add_cell_data(cell_data, split_pattern=split_pattern)

        # Save annotated object
        output_path = os.path.join(output_dir, f"{sample_id}_cistopic_obj.pkl")
        with open(output_path, "wb") as f:
            pickle.dump(cistopic_obj, f)
        print(f"Saved annotated CistopicObject: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add cell metadata to CistopicObjects")
    parser.add_argument("--cistopic_objects_pickle", required=True,
                        help="Pickle file containing list of CistopicObjects")
    parser.add_argument("--cell_data", required=True,
                        help="TSV file with cell metadata (barcodes as rows)")
    parser.add_argument("--output_dir", required=True,
                        help="Directory to save annotated CistopicObjects")
    parser.add_argument("--split_pattern", default="-",
                        help="Pattern used to split barcodes (default '-')")

    args = parser.parse_args()
    main(args.cistopic_objects_pickle, args.cell_data, args.output_dir, args.split_pattern)

