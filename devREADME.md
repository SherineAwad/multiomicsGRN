# Clustering and Visualization Step in scRNA-seq

## Overview

This step performs **cell-level preprocessing, dimensionality reduction, clustering, and visualization** for single-cell RNA-seq data.  

The goal is to **group similar cells together** based on their gene expression profiles and prepare the data for downstream analyses and interpretation.

**Key conceptual operations:**

1. **Normalization & log transformation**  
   - Makes counts comparable across cells and reduces the influence of highly variable cells.

2. **Dimensionality reduction (PCA)**  
   - Condenses thousands of genes into a smaller set of principal components capturing most variation.

3. **Neighbor graph construction**  
   - Identifies cells with similar expression profiles in the reduced-dimensional space.

4. **Clustering (Leiden algorithm)**  
   - Groups similar cells together, each cluster potentially representing a distinct cell type or state.

5. **Visualization (UMAP)**  
   - Projects cells into 2D space to see cluster structure and relationships.

6. **Marker gene plotting**  
   - Evaluates the activity of known marker genes to interpret cluster identities.

---

## Inputs

1. **AnnData object (`.h5ad`)**  
   - Contains the single-cell gene expression matrix.

2. **Marker gene list**  
   - Text file of known genes to evaluate cluster identity.

---

## Outputs

1. **Clustered AnnData object**  
   - Includes UMAP coordinates, cluster labels (`leiden`), and sample metadata.

2. **UMAP plots**  
   - Cells colored by clusters and by sample.

3. **Dotplots of marker genes**  
   - Summarizes marker expression across clusters.

4. **UMAPs for individual marker genes**  
   - Visualizes expression of each marker gene across cells.

---

## Why This Step is Recommended

1. **Defines cell types or states**  
   - Clustering organizes cells based on similarity, providing a clear structure for interpretation.

2. **Supports downstream analyses**  
   - Clusters and cell type labels can guide pseudobulk aggregation, marker gene evaluation, or further integrative analyses.

3. **Improves interpretability**  
   - Visualization and marker gene plotting help assign biological meaning to clusters, making results easier to understand and communicate.



# Pycistopic 


## Preprocessing Step in scATAC-seq 

## Overview

This step is the **“pseudobulk preprocessing”** for scATAC-seq data.
Single-cell ATAC-seq is very sparse — most genomic regions have zero reads in individual cells. Analyzing each cell separately for peak calling is noisy and computationally intensive.

**Pseudobulk preprocessing does the following:**

- Aggregates reads from multiple cells of the same type into a single profile.
- Generates clean, genome-wide coverage tracks for each cell type and sample.
- Creates a “virtual bulk ATAC-seq experiment” for each group of cells.

---

## Inputs

1. **Fragment files**: Raw scATAC-seq fragments per sample, containing coordinates of all sequenced fragments per cell.
2. **Cell metadata**: Maps each barcode to a cell type and sample.
3. **Chromosome sizes**: Defines the genome boundaries for coverage tracks.

---

## Outputs

1. **BED files**: For each cell type × sample, showing all fragment coordinates. These are used later for peak calling.
2. **BigWig files**: Normalized genome coverage tracks, useful for visualization and downstream analyses.
3. **TSV reference files**: Lists of the generated BED and BigWig files for easy tracking.

---

## Analogy

Imagine each cell as a noisy sensor detecting “open windows” (accessible regions) in a city (the genome).
Each sensor alone is unreliable, but if you **combine all sensors of the same type**, you get a clear map of which windows are open for that cell type. That combined map is your **pseudobulk**.


# Pseudobulk Export Step in pycisTopic

## Overview

This step generates **aggregated pseudobulk profiles** from single-cell ATAC-seq data to create clean coverage tracks for downstream analysis.  

- Single-cell ATAC data is sparse, so combining fragments from cells of the same **cell type × sample** produces more robust profiles.  
- This step **takes the BED and BigWig files generated during preprocessing** as input.  

---

## Inputs

1. **BED and BigWig files from preprocessing**  
   - Represent aggregated fragment coordinates and coverage per cell type × sample.  

2. **Barcode-to-celltype metadata**  
   - Maps each cell to its cell type and sample, ensuring correct assignment for aggregation.  

3. **Chromosome sizes**  
   - Defines genome boundaries for coverage tracks.

---

## What it does?

1. **Map fragments to cell types**  
   - Assigns fragments from the BED/BigWig files to the correct cell type × sample.  

2. **Aggregate fragments**  
   - Combines all fragments for each group into a single pseudobulk profile.  

3. **Generate updated coverage tracks**  
   - **BED files:** lists of fragments per pseudobulk.  
   - **BigWig files:** normalized coverage for visualization and downstream analysis.  

4. **Save reference files**  
   - Stores paths to all generated BED and BigWig files for later use.

---

## Outputs

1. **BED files** – aggregated fragment coordinates per cell type × sample.  
2. **BigWig files** – normalized genome coverage tracks per pseudobulk.  
3. **TSV reference files** – lists of all BED and BigWig file paths.

---

## Analogy

Each cell is a noisy “sensor” detecting open windows in the genome.  
- This step **combines all sensors of the same type into a single, clear map**, producing robust accessibility profiles for each cell type × sample.




# Peak Calling Step with MACS2

## Overview

This step identifies **peaks**, i.e., genomic regions that are significantly enriched for ATAC-seq fragments.  

- Pseudobulk files summarize the signal for each **cell type × sample**.  
- MACS2 uses these aggregated fragment coordinates to detect regions with **high accessibility**, likely representing regulatory elements such as enhancers or promoters.  
- The output peaks are later used to generate **consensus peaks** and for topic modeling in pycisTopic.

---

## Inputs

1. **Pseudobulk BED files** from the previous step  
   - Contain aggregated fragment coordinates per cell type × sample.

2. **Genome size** (e.g., `mm` for mouse)  
   - Required for MACS2 to estimate background signal and calculate statistical significance.

---

## What it does?: 

1. **Convert fragment files to BED format**  
   - Ensures compatibility with MACS2 (chromosome, start, end).  

2. **Call peaks using MACS2**  
   - Scans the genome to find regions with **more fragments than expected by chance**.  
   - Parameters like `--shift` and `--extsize` adjust fragment positioning to capture accessible regions accurately.  
   - Each pseudobulk sample is processed independently.  

3. **Store peak files**  
   - Produces a BED-like file of called peaks for each pseudobulk.

---

## Outputs

1. **Peak files (`.narrowPeak`)**  
   - Lists genomic regions with significantly enriched ATAC fragments per pseudobulk.  

2. **BED files of peaks (optional)**  
   - Can be used for merging into consensus peaks across samples.  

3. **Auxiliary MACS2 files**  
   - Statistics and logs for quality control (e.g., peak scores, fragment pileups).

---

## Analogy

Think of the pseudobulk BED as a **heat map of open windows across the city**.  
- MACS2 finds **clusters of “hot spots”** where many windows are open at once — these are your peaks.  
- Each peak represents a genomic region with strong evidence of accessibility in that cell type.

---

✅ **Summary**

- **Purpose:** Detect high-confidence accessible regions from aggregated pseudobulk data.  
- **Inputs:** Pseudobulk BED files, genome size.  
- **Outputs:** MACS2 peak files (BED/narrowPeak) per cell type × sample, ready for **consensus peak generation** and downstream analysis.


# Consensus Peak Generation Step

## Overview

This step creates a **consensus peak set** by merging peaks called by MACS2 across all pseudobulk samples.  

- Each pseudobulk sample has its own set of peaks (from MACS2).  
- To perform downstream analyses (e.g., cistopic topic modeling, DAR analysis), you need a **single unified set of peaks** representing all accessible regions across samples.  

---

## Inputs

1. **MACS2 peak files (`.narrowPeak`)**  
   - Generated in the previous MACS2 step for each pseudobulk sample.  

2. **Combined BED file (optional intermediate)**  
   - A concatenation of all individual narrowPeak files before merging.  

---

# Merge  

1. **Combine all MACS2 peak files**  
   - Concatenate narrowPeak files from all pseudobulk samples into a single file.  

2. **Sort and merge overlapping peaks**  
   - Using `bedtools sort` and `bedtools merge`, overlapping or adjacent peaks are merged into single regions.  
   - This produces a **non-redundant consensus peak set**, capturing all regions identified in any sample.  

---

## Outputs

1. **Consensus peaks BED file**  
   - A single BED file containing all merged peaks across samples.  
   - This is used for creating **cistopic objects** and for downstream analyses like differential accessibility.  

---

## Analogy

Imagine each pseudobulk sample as a city map showing open windows (peaks).  
- Each map is slightly different because some windows are only open in certain samples.  
- Merging all maps produces a **master map of all windows ever open** — this is the consensus peak set.

---

✅ **Summary**

- **Purpose:** Generate a unified, non-redundant set of peaks across all pseudobulk samples.  
- **Inputs:** MACS2 narrowPeak files from pseudobulk samples.  
- **Outputs:** Consensus peak BED file, ready for cistopic object creation and downstream analyses.




