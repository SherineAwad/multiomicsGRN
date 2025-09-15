#!/usr/bin/env python3
import os
import argparse
import pycisTopic as pct

def main(genome_name: str, output_file: str):
    """
    Generate a TSS BED file using pycisTopic gene annotations.
    """
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Generate TSS annotation from Ensembl
    print(f"[INFO] Fetching TSS annotations for {genome_name}...")
    tss_annotation_df = pct.gene_annotation.get_tss_annotation_from_ensembl(
        biomart_name=genome_name
    )

    # Write TSS annotation to BED file
    pct.gene_annotation.write_tss_annotation_to_bed(
        tss_annotation_bed_df_pl=tss_annotation_df,
        tss_annotation_bed_filename=output_file
    )

    print(f"[INFO] TSS BED file saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate TSS BED file for a given genome using pycisTopic")
    parser.add_argument(
        "-g", "--genome", 
        default="mmusculus_gene_ensembl", 
        help="Ensembl genome name (default: mmusculus_gene_ensembl)"
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Output BED file path"
    )
    args = parser.parse_args()
    main(args.genome, args.output)

