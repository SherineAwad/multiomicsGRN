#!/bin/bash
# File: make_genome_annotation.sh

GTF="/nfs/turbo/umms-thahoang/sherine/REFERENCES/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf"
OUT="genome_annotation.tsv"

# Remove header lines starting with #
# Extract either "gene" or "transcript" features
# Columns: chrom  start  end  strand  gene_id  gene_name
awk 'BEGIN{OFS="\t"} 
     $0 !~ /^#/ && ($3=="gene" || $3=="transcript") {
         match($0, /gene_id "([^"]+)"/, gid);
         match($0, /gene_name "([^"]+)"/, gname);
         if (gid[1]!="" && gname[1]!="") print $1, $4, $5, $7, gid[1], gname[1]
     }' "$GTF" > "$OUT"

echo "SCENIC+-friendly genome_annotation.tsv created at $OUT"

