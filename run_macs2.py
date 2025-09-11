#!/usr/bin/env python3
import os
import subprocess
import glob

# Input and output directories
frag_dir = "./consensus_peak_calling/pseudobulk_bed_files"
out_dir = "./consensus_peak_calling/macs2_peaks"
genome = "mm"  # mouse, change to "hs" for human

os.makedirs(out_dir, exist_ok=True)

# Find all fragment files
frag_files = glob.glob(os.path.join(frag_dir, "*.fragments.tsv.gz"))

for frag in frag_files:
    sample = os.path.basename(frag).replace(".fragments.tsv.gz", "")
    print(f">>> Processing {sample}")

    bed_file = os.path.join(out_dir, f"{sample}.bed")

    # Convert fragments to BED
    cmd_bed = f"zcat {frag} | awk 'BEGIN{{OFS=\"\\t\"}} {{print $1,$2,$3}}' | sort -k1,1 -k2,2n > {bed_file}"
    subprocess.run(cmd_bed, shell=True, check=True)

    # Call MACS2
    cmd_macs = [
        "macs2", "callpeak",
        "-t", bed_file,
        "-f", "BED",
        "-n", sample,
        "--outdir", out_dir,
        "--nomodel",
        "--shift", "-100",
        "--extsize", "200",
        "-g", genome
    ]
    subprocess.run(cmd_macs, check=True)

    print(f">>> Done {sample}")

