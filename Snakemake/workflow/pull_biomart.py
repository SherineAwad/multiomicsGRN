#!/usr/bin/env python3
import pandas as pd
import requests
import gzip
import io
import subprocess
import sys

def download_ucsc_genes():
    """Download gene annotations from UCSC in SCENIC+ format"""
    print("Downloading gene annotations from UCSC...")
    
    # Download UCSC genes
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz"
    response = requests.get(url)
    
    with gzip.open(io.BytesIO(response.content), 'rt') as f:
        lines = f.readlines()
    
    # Parse GTF and convert to SCENIC+ format with ALL required columns
    genes = []
    for line in lines:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if parts[2] == 'transcript':
            chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            
            # Parse attributes to get gene name
            attrs = {}
            for item in parts[8].split(';'):
                if ' ' in item.strip():
                    key, value = item.strip().split(' ', 1)
                    attrs[key] = value.strip('"')
            
            gene_name = attrs.get('gene_name', '') or attrs.get('gene_id', '')
            
            # Calculate Transcription Start Site based on strand
            if strand == '+':
                tss = start - 1  # 0-based, TSS is start for + strand
            else:
                tss = end  # TSS is end for - strand
            
            # Use ALL column names SCENIC+ expects
            genes.append({
                'Chromosome': chrom,
                'Start': start - 1,  # Convert to 0-based
                'End': end,
                'Strand': strand,
                'Gene': gene_name,
                'Transcription_Start_Site': tss
            })
    
    df = pd.DataFrame(genes)
    # Ensure we have ALL the exact columns SCENIC+ wants
    df = df[['Chromosome', 'Start', 'End', 'Strand', 'Gene', 'Transcription_Start_Site']]
    df.to_csv('genome_annotation.tsv', sep='\t', index=False)
    print(f"Saved {len(df)} genes to genome_annotation.tsv")
    print("Columns:", df.columns.tolist())

def download_chromsizes():
    """Download chromosome sizes and convert to SCENIC+ format"""
    print("Downloading chromosome sizes from UCSC...")
    
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
    df = pd.read_csv(url, sep='\t', header=None, names=['Chromosome', 'Size'])
    
    # Convert to SCENIC+ format: Chromosome, Start, End
    # For chromsizes, Start is always 0, End is the chromosome size
    df['Start'] = 0
    df['End'] = df['Size']
    
    # Keep only the columns SCENIC+ expects
    df = df[['Chromosome', 'Start', 'End']]
    df.to_csv('chromsizes.tsv', sep='\t', index=False)
    print(f"Saved {len(df)} chromosomes to chromsizes.tsv")
    print("Columns:", df.columns.tolist())

def main():
    """Download and format all files for SCENIC+"""
    try:
        download_ucsc_genes()
        download_chromsizes()
        print("✅ SUCCESS: All files downloaded and formatted for SCENIC+")
        
        # Verify files
        print("\n=== Verification ===")
        print("genome_annotation.tsv:")
        subprocess.run(['head', '-n', '3', 'genome_annotation.tsv'])
        print("\nchromsizes.tsv:")
        subprocess.run(['head', '-n', '3', 'chromsizes.tsv'])
        
    except Exception as e:
        print(f"❌ ERROR: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
