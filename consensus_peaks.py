#!/usr/bin/env python3
import os
import glob
import subprocess
import argparse

# -----------------------------
# 1. Define directories
# -----------------------------
macs2_dir = "./consensus_peak_calling/macs2_peaks"
consensus_out = "./consensus_peak_calling/consensus"
os.makedirs(consensus_out, exist_ok=True)
parser = argparse.ArgumentParser()
parser.add_argument("combined_bed", help="Path to combined bed file")
parser.add_argument("consensus_peaks", help="Path to consensus peaks bed file")
args = parser.parse_args()
combined_bed = args.combined_bed
consensus_peaks = args.consensus_peaks 
# -----------------------------
# 2. Merge all narrowPeak files
# -----------------------------
# Get all MACS2 narrowPeak files
narrow_files = glob.glob(os.path.join(macs2_dir, "*.narrowPeak"))

# Combine all narrowPeak files into one
combined_narrow = os.path.join(consensus_out, combined_bed)
with open(combined_narrow, "w") as out_f:
    for f in narrow_files:
        with open(f) as in_f:
            for line in in_f:
                out_f.write(line)

# Sort and merge peaks using bedtools
consensus_bed = os.path.join(consensus_out, consensus_peaks)
subprocess.run(f"bedtools sort -i {combined_narrow} | bedtools merge > {consensus_bed}",
               shell=True, check=True)
print(f"Consensus peak set created: {consensus_bed}")


