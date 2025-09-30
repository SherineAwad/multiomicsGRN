# 🔴🔴🔴 PART A:  Preprcessing of scRNA -  Clustering and Visualization

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

## Part A Results: scRNA-seq results    

![Before Filtering QC metrics](figures/violin_QC.png)
Violin plots displaying quality control metrics such as number of genes detected per cell, total counts, and percentage of mitochondrial gene expression.

### Filtering Criteria

Quality control filtering was applied to exclude low-quality cells and potential doublets. Cells were retained only if they met all the following conditions:

- Number of genes detected per cell between **500 and 7000**
- Total counts per cell between **1000 and 30,000**
- Percentage of mitochondrial gene counts less than **25%**

This filtering step ensures removal of dead or dying cells, doublets, and technical artifacts, thereby improving the quality of downstream analyses.


### Additional Analysis Figure
![After Filtering QC metrics](figures/violin_AfterQC.png)

### UMAP and Clusters UMAP


![Samples UMAP](figures/umap_clustered_mNeurog2_Samples.png)


![Clusters UMAP](figures/umap_clustered_mNeurog2_Clusters.png)

### DOTPLOT 
![Dotplot](figures/clustered_mNeurog2_Dotplot.png)

### Marker Genes UMAP

<img src="figures/umap_clustered_mNeurog2_Abca8a.png?v=4" alt="Abca8a" width="33%"><img src="figures/umap_clustered_mNeurog2_Ccr2.png?v=4" alt="Ccr2" width="33%"><img src="figures/umap_clustered_mNeurog2_Hes1.png?v=4" alt="Hes1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Notch1.png?v=4" alt="Notch1" width="33%"><img src="figures/umap_clustered_mNeurog2_Rlbp1.png?v=4" alt="Rlbp1" width="33%"><img src="figures/umap_clustered_mNeurog2_Acta2.png?v=4" alt="Acta2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Chat.png?v=4" alt="Chat" width="33%"><img src="figures/umap_clustered_mNeurog2_Hes5.png?v=4" alt="Hes5" width="33%"><img src="figures/umap_clustered_mNeurog2_Nrl.png?v=4" alt="Nrl" width="33%">
<img src="figures/umap_clustered_mNeurog2_Rpe65.png?v=4" alt="Rpe65" width="33%"><img src="figures/umap_clustered_mNeurog2_Apoe.png?v=4" alt="Apoe" width="33%">
<img src="figures/umap_clustered_mNeurog2_Insm1.png?v=4" alt="Insm1" width="33%"><img src="figures/umap_clustered_mNeurog2_Olig2.png?v=4" alt="Olig2" width="33%"><img src="figures/umap_clustered_mNeurog2_Sebox.png?v=4" alt="Sebox" width="33%">
<img src="figures/umap_clustered_mNeurog2_Aqp4.png?v=4" alt="Aqp4" width="33%"><img src="figures/umap_clustered_mNeurog2_Csf1r.png?v=4" alt="Csf1r" width="33%"><img src="figures/umap_clustered_mNeurog2_Isl1.png?v=4" alt="Isl1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Otx2.png?v=4" alt="Otx2" width="33%"><img src="figures/umap_clustered_mNeurog2_Slc17a7.png?v=4" alt="Slc17a7" width="33%"><img src="figures/umap_clustered_mNeurog2_Arr3.png?v=4" alt="Arr3" width="33%">
<img src="figures/umap_clustered_mNeurog2_Elavl3.png?v=4" alt="Elavl3" width="33%"><img src="figures/umap_clustered_mNeurog2_Kcnj8.png?v=4" alt="Kcnj8" width="33%"><img src="figures/umap_clustered_mNeurog2_Pax2.png?v=4" alt="Pax2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Slc1a3.png?v=4" alt="Slc1a3" width="33%"><img src="figures/umap_clustered_mNeurog2_Ascl1.png?v=4" alt="Ascl1" width="33%"><img src="figures/umap_clustered_mNeurog2_Elavl4.png?v=4" alt="Elavl4" width="33%">
<img src="figures/umap_clustered_mNeurog2_Lhx1.png?v=4" alt="Lhx1" width="33%"><img src="figures/umap_clustered_mNeurog2_Pax6.png?v=4" alt="Pax6" width="33%"><img src="figures/umap_clustered_mNeurog2_Slc6a9.png?v=4" alt="Slc6a9" width="33%">
<img src="figures/umap_clustered_mNeurog2_Atoh7.png?v=4" alt="Atoh7" width="33%"><img src="figures/umap_clustered_mNeurog2_Emx1.png?v=4" alt="Emx1" width="33%"><img src="figures/umap_clustered_mNeurog2_Lhx2.png?v=4" alt="Lhx2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Pou4f2.png?v=4" alt="Pou4f2" width="33%"><img src="figures/umap_clustered_mNeurog2_Sox11.png?v=4" alt="Sox11" width="33%"><img src="figures/umap_clustered_mNeurog2_Bsn.png?v=4" alt="Bsn" width="33%">
<img src="figures/umap_clustered_mNeurog2_Foxn4.png?v=4" alt="Foxn4" width="33%"><img src="figures/umap_clustered_mNeurog2_Lhx4.png?v=4" alt="Lhx4" width="33%"><img src="figures/umap_clustered_mNeurog2_Prdm1.png?v=4" alt="Prdm1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Sox9.png?v=4" alt="Sox9" width="33%"><img src="figures/umap_clustered_mNeurog2_Cabp5.png?v=4" alt="Cabp5" width="33%"><img src="figures/umap_clustered_mNeurog2_Gad1.png?v=4" alt="Gad1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Malat1.png?v=4" alt="Malat1" width="33%"><img src="figures/umap_clustered_mNeurog2_Prdx6.png?v=4" alt="Prdx6" width="33%"><img src="figures/umap_clustered_mNeurog2_Tfap2a.png?v=4" alt="Tfap2a" width="33%">
<img src="figures/umap_clustered_mNeurog2_Calb1.png?v=4" alt="Calb1" width="33%"><img src="figures/umap_clustered_mNeurog2_Gfap.png?v=4" alt="Gfap" width="33%"><img src="figures/umap_clustered_mNeurog2_mt-Atp6.png?v=4" alt="mt-Atp6" width="33%">
<img src="figures/umap_clustered_mNeurog2_Rbfox3.png?v=4" alt="Rbfox3" width="33%"><img src="figures/umap_clustered_mNeurog2_Tie1.png?v=4" alt="Tie1" width="33%"><img src="figures/umap_clustered_mNeurog2_Calb2.png?v=4" alt="Calb2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Glul.png?v=4" alt="Glul" width="33%"><img src="figures/umap_clustered_mNeurog2_Neurog2.png?v=4" alt="Neurog2" width="33%"><img src="figures/umap_clustered_mNeurog2_Rho.png?v=4" alt="Rho" width="33%">
<img src="figures/umap_clustered_mNeurog2_Vim.png?v=4" alt="Vim" width="33%">

### Annotations 

![Annotations ON](figures/annotated_clustered_mNeurog2_annotationsON.png)

![Annotations](figures/annotated_clustered_mNeurog2_annotations.png)

### Cell counts 

This is the cell counts in clusters based on barcodes:

| Cluster      | Cell Count |
|--------------|-----------|
| MG           | 3206      |
| BC           | 1298      |
| Rod          | 997       |
| MGPC         | 612       |
| Microglia    | 344       |
| AC           | 103       |
| Cone         | 48        |
| Astrocyte    | 4         |
|**Total cells** | **6612** |


# 🔴🔴🔴 PART B: Pycistopic: preprocessing the ATAC part and integrating metadata from scRNA 

# 1. Preprocessing scATAC-seq 

## Overview

This step frames the **challenge of scATAC-seq data** and introduces the solution of pseudobulk aggregation.  
- Single-cell ATAC-seq is **extremely sparse** — most individual cells have very few reads across the genome.  
- Calling peaks or interpreting accessibility patterns on a cell-by-cell basis is unreliable.  
- To overcome this, we **aggregate reads across cells of the same type and condition**, creating **pseudobulk profiles**.  

This gives us a foundation for robust downstream analysis, much like turning many weak sensors into one strong signal.

---

## What does this achieve?

- **Stabilizes the signal**: reduces noise by pooling reads from similar cells.  
- **Enables peak calling**: MACS2 and other tools require sufficient depth, which pseudobulk provides.  
- **Improves interpretability**: each pseudobulk resembles a bulk ATAC-seq experiment, making results biologically clearer.  

---

## Inputs 

1. **Fragment files**  
   - Raw outputs from the scATAC-seq experiment (e.g., Cell Ranger ATAC).  

2. **Cell metadata**  
   - Barcode annotations (e.g., cell type, sample) generated upstream using Seurat, Scanpy, or similar tools.  

3. **Chromosome sizes**  
   - Genome reference information needed for building coverage tracks.  

---

## Outputs 

- **Aggregated accessibility profiles** per cell type × sample are **saved as TSV tables**.  
  - Each TSV lists the fragment coordinates for the pseudobulk group.  
  - These TSV files are **intermediate outputs**, not yet in BED/BigWig format.  
  - The TSVs **define what will later be converted into BED and BigWig files** for downstream analysis.

---

## Analogy

Imagine many faint, crackling radios (single cells) all trying to broadcast the same station.  
Individually, the signal is unclear — but by combining them, you reconstruct the song clearly.  
The **TSV tables** are like the combined raw audio recordings — you now have the data ready to be processed into a clean, playable track.  
i

## scATAC preanalyis by Seurat as input to pycistopic preprocessing 

![ATAC UMAP](ATAC_samples.png)

### Clustering 

![ATAC CLUSTERS](ATAC_clusters.png)


---

# 2. Pseudobulk Export with pycisTopic 

## Overview

This step **executes the pseudobulk aggregation** using `pycisTopic`.  
- Takes the conceptual idea and TSV intermediates from Step 1 and **produces the actual BED/BigWig files** for downstream analysis.  
- Aggregates reads per **cell type × sample**, generating clean coverage tracks for visualization and peak calling.  

---

## What does it do?

1. **Maps fragments to annotated cells**  
   - Uses metadata (from Seurat/Scanpy) to assign each fragment to its correct cell type × sample group.  

2. **Aggregates fragments into pseudobulk TSVs** (if not already created)  
   - Combines all reads for each group into a single table of fragment coordinates.  

3. **Converts TSVs to standard formats**  
   - **BED files:** fragment positions formatted for downstream peak calling (MACS2, etc.).  
   - **BigWig files:** normalized coverage tracks for visualization in genome browsers.  

4. **Organizes outputs**  
   - Saves **TSV reference files** with paths to all BED and BigWig files for downstream pipelines.  

---

## Inputs 

1. **Fragment files**: `fragments.tsv.gz` per sample (from the experiment).  
2. **Cell metadata CSV**: barcode → cell type/sample mapping (from Seurat/Scanpy).  
3. **Chromosome sizes**: defines genome boundaries (from UCSC or Ensembl).  

---

## Outputs

1. **BED files** – fragment coordinates per pseudobulk (cell type × sample).  
2. **BigWig files** – normalized coverage tracks for genome browsers or QC.  
3. **TSV lists** – references of all generated BED/BigWig paths.  

---

## Analogy

Step 1 was like **planning and creating the combined raw audio (TSV tables)**.  
Step 2 is where we **mix, normalize, and export it into clean, usable formats (BED/BigWig)** for analysis.  


# 3. Peak Calling Step with MACS2

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

## What it does? 

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

After pseudobulk aggregation and MACS2 peak calling, each sample (or cell type × sample) has its **own set of peaks**.  
- Different samples may have slightly different peaks because of biological variability or sequencing depth.  
- To do comparative analyses across samples or feed data into pycisTopic, we need **one unified set of peaks** — the **consensus peak set**.  

Think of it as **finding all regions that are open in at least one sample and merging overlapping regions into a master list**.

---

## What does it do?

1. **Collect all individual MACS2 peaks**  
   - Each pseudobulk BED file from Step 2 has peaks called separately using MACS2.  

2. **Combine all peak files**  
   - Concatenate all peak files into one large table.  

3. **Sort and merge overlapping peaks**  
   - Using tools like **bedtools sort** and **bedtools merge**.  
   - If two peaks from different samples overlap, they are merged into a single peak region.  
   - The result is a **non-redundant, genome-wide consensus peak set**.

---

## Inputs

1. **MACS2 peak files** (`*.narrowPeak` or BED)  
   - Generated per sample/pseudobulk in the previous MACS2 step.  

2. **Optional parameters** (if any) for merging, e.g., minimum overlap or padding.  

---

## Outputs

1. **Consensus peaks BED file**  
   - One unified list of all peak regions across samples.  
   - This file is later used to define **features/regions in pycisTopic objects**.  
2. **Intermediate combined BED** (optional)  
   - The unsorted concatenated peaks before merging.  

---

## Why it is needed

- Ensures **all samples/cell types are analyzed on the same set of genomic regions**, which is critical for:  
  - Accurate topic modeling (pycisTopic LDA)  
  - Differential accessibility analysis  
  - Integration with scRNA-seq metadata  

---

## Analogy

Imagine several people drawing maps of the same city, each highlighting where the windows are open.  
- Each map may differ slightly.  
- The consensus peak step **merges all maps into one master map** that shows all open windows observed across everyone.  


---

# 5. Merge  

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



# 6. TSS Generation Step in pycisTopic

## Overview

This step generates a **BED file containing the transcription start sites (TSSs)** of genes for the reference genome.  

- TSS regions are important for **quality control**, such as checking **TSS enrichment** in ATAC-seq data.  
- It can also be used to **annotate peaks** with nearby genes for downstream analyses.  

---

## Inputs

1. **Reference genome annotation**  
   - Provided by pycisTopic via Ensembl gene annotations (e.g., `"mmusculus_gene_ensembl"`).  
   - Contains gene coordinates, including TSS positions.

2. **Genome build specification**  
   - UCSC or Ensembl genome coordinates (e.g., `"mm10"` for mouse).  

---

## What it does?

1. **Extract TSS positions from gene annotations**  
   - For each gene, determine the **chromosome, strand, and start position** corresponding to the TSS.  

2. **Convert to BED format**  
   - Each TSS is represented as a genomic interval suitable for downstream QC or annotation.

3. **Save TSS BED file**  
   - The BED file can be used for **TSS enrichment analysis** or annotating peaks to genes.

---

## Outputs

1. **TSS BED file**  
   - Contains genomic coordinates of all transcription start sites for the reference genome.  
   - Example path: `outs/qc/tss_mm10.bed`.

---

## Connection to Previous Steps

- This step **does not directly depend on the pseudobulk or MACS2 peaks**.  
- It **uses the reference genome annotation**, so it is independent of previous sample-specific steps.  
- The generated TSS BED file is primarily used for **quality control** (e.g., checking enrichment of ATAC signal at TSSs) and for **peak annotation downstream**.

---

✅ **Summary**

- **Purpose:** Generate a reference set of TSSs in BED format.  
- **Inputs:** Reference genome annotation, genome build.  
- **Outputs:** TSS BED file, used for QC and peak annotation.  
- **Relation to previous steps:** Independent; serves as a reference for evaluating and annotating pseudobulk/MACS2 peaks.


# 7. Creating Cistopic Objects Step

## Overview

This step creates a **cistopic object**, which is the central data structure used by pycisTopic for **topic modeling of chromatin accessibility**.  

- The cistopic object organizes **fragment data, peak regions, and QC information** in a way suitable for downstream analyses.  
- It is essentially a **single-cell peak-by-cell matrix** stored in a Python pickle file, with metadata attached.  

---

## Inputs

1. **Fragments dictionary (`fragments_dict.pkl`)**  
   - Maps each sample name to its fragment file (`.tsv.gz`).  
   - Produced or referenced from earlier steps (preprocess/pseudobulk fragment info).  

2. **QC results (`qc_barcodes_thresholds.pkl`)**  
   - Information about filtered/high-quality cells, such as thresholds for minimum fragments, TSS enrichment, etc.  
   - Derived from QC steps performed earlier.  

3. **Consensus peaks BED file**  
   - Generated in the **previous consensus peak step**.  
   - Defines the genomic regions (peaks) that will form the rows of the cell-by-peak matrix.  

4. **Blacklist BED file (`mm10-blacklist.v2.bed`)**  
   - Contains regions of the genome known to produce artefactual signals, which are excluded.  

5. **QC output directory**  
   - Provides additional QC-related files for integration.  

6. **Number of CPUs (`n_cpu`)**  
   - Used for parallel processing.  

---

## What it does?

1. **Load fragments and QC metadata**  
   - Filter cells based on QC thresholds.  

2. **Map fragments to consensus peaks**  
   - Assign each fragment from a cell to a peak, creating a **peak-by-cell matrix**.  

3. **Remove blacklisted regions**  
   - Exclude problematic regions that may generate false positives.  

4. **Store in cistopic object**  
   - Includes:
     - Cell-by-peak matrix  
     - Cell metadata (sample, cell type, QC info)  
     - Peak metadata (genomic coordinates, annotation)  

5. **Save as Python pickle**  
   - This object will be used for **topic modeling, clustering, and downstream analysis**.

---

## Outputs

1. **Cistopic object (`cistopic_objects_mm10.pkl`)**  
   - Single Python object storing the filtered peak-by-cell matrix and associated metadata.  
   - Ready for running **topic modeling (LDA)** and other analyses in pycisTopic.

---

## Connection to Previous Steps

- **Takes consensus peaks** (from merged MACS2 outputs).  
- **Uses fragment files** (from preprocessing/pseudobulk).  
- **Incorporates QC information** (filtered cells, TSS enrichment, etc.).  

- This step is **critical** because it transforms raw and pseudobulk fragment data into a structured object suitable for all downstream pycisTopic analyses.

## QC plots of our cistopic object 


![General QC TH1](outs/qc_plots/TH1_qc.png)

![General QC TH2](outs/qc_plots/TH2_qc.png)

![Barcode QC TH1](outs/barcode_qc_plots/TH1_barcode_qc.png)

![Barcode QC TH2](outs/barcode_qc_plots/TH2_barcode_qc.png)


# 8. Merging Cistopic Objects Step

## Overview

This step merges **one or more cistopic objects** into a single unified object.  

- In workflows with multiple samples, batches, or preprocessing runs, each cistopic object may represent a separate sample or subset of cells.  
- Merging combines them into a **single cistopic object**, making downstream analyses (topic modeling, clustering, DAR analysis) easier and consistent across all cells.  

---

## Inputs

1. **Cistopic object(s) (`cistopic_objects_mm10.pkl`)**  
   - Created in the previous step.  
   - Contains:
     - Peak-by-cell matrices  
     - Cell metadata (cell type, sample, QC info)  
     - Peak metadata  

2. **Output path (`merged_cistopic.pkl`)**  
   - Where the merged cistopic object will be saved.  

---

## What it does?

1. **Load all input cistopic objects**  
   - One object per sample or batch.  

2. **Combine peak matrices**  
   - Align peaks across objects to ensure **consistent genomic coordinates**.  

3. **Merge cell metadata**  
   - Keep track of **cell type, sample, QC info** for all cells.  

4. **Create unified cistopic object**  
   - Contains all cells and peaks in a single structure.  

5. **Save merged object**  
   - Stored as a Python pickle for use in downstream pycisTopic steps.

---

## Outputs

1. **Merged cistopic object (`merged_cistopic.pkl`)**  
   - Single object containing all cells and peaks.  
   - Ready for:
     - LDA/topic modeling  
     - Clustering  
     - Differential accessibility analysis  

---

## Connection to Previous Steps

- **Takes individual cistopic objects** generated in the `create_cistopic_objects.py` step.  
- Aligns and merges them so that **all samples/cells are in one unified object** for downstream pycisTopic analysis.  

✅ **Summary:**  
- **Purpose:** Create a single unified cistopic object from multiple inputs.  
- **Inputs:** One or more cistopic objects (pickle files).  
- **Outputs:** Merged cistopic object ready for topic modeling and downstream analyses.



# 9. Adding scRNA-seq Metadata to Cistopic Objects

## Overview

This step integrates **scRNA-seq-derived metadata** into the merged cistopic object.  

- Single-cell RNA-seq preprocessing (clustering, cell type annotation) provides **cell type labels, sample IDs, or other metadata**.  
- Attaching this information to the cistopic object allows **linking chromatin accessibility topics to known cell types** for interpretation.  

---

## Inputs

1. **Merged cistopic object (`merged_cistopic.pkl`)**  
   - Created in the previous merge step.  
   - Contains the **ATAC-seq peak-by-cell matrix** and ATAC metadata.  

2. **scRNA-seq metadata CSV (`scRNA_barcodes.csv`)**  
   - Contains cell barcodes and associated information such as:
     - Cell type labels  
     - Sample ID  
     - Cluster assignments  

3. **Output pickle path (`merged_with_meta.pkl`)**  
   - Where the updated cistopic object with metadata will be saved.  

---

## What it does?

1. **Load merged cistopic object**  
   - Contains all ATAC-seq cells in a single structure.  

2. **Load scRNA metadata**  
   - Map cell barcodes to those in the cistopic object.  

3. **Attach metadata to cistopic object**  
   - Adds new columns in the cell metadata (obs) of the cistopic object.  
   - This allows **annotating cells with known cell types or clusters** derived from scRNA-seq.  

4. **Save updated cistopic object**  
   - The new object contains both chromatin accessibility data and scRNA-derived annotations.  

---

## Outputs

1. **Cistopic object with scRNA metadata (`merged_with_meta.pkl`)**  
   - Each cell now has ATAC-seq data **and** associated scRNA-seq-derived labels.  
   - Enables cell type–specific analyses in topic modeling, DAR identification, and visualization.  

---

## Connection to Previous Steps

- **Takes the merged cistopic object** from the previous merging step.  
- **Takes metadata from scRNA preprocessing**, including clusters and cell type labels.  
- **Purpose:** Link ATAC-seq profiles to known cell types for biological interpretation and downstream analyses.  

✅ **Summary:**  
- **Why needed:** Without scRNA metadata, ATAC-seq topics cannot be directly assigned to cell types.  
- **Input:** Merged cistopic object + scRNA metadata.  
- **Output:** Annotated cistopic object ready for topic modeling and cell type–specific analysis.


# 10. Topic Modeling with Mallet (run_mallet.py)

## Overview

This step performs **Latent Dirichlet Allocation (LDA) topic modeling** on the chromatin accessibility data stored in the cistopic object.  

- Each **topic represents a set of genomic regions (peaks) that tend to be accessible together** across cells.  
- Topic modeling reduces the high-dimensional peak-by-cell matrix into a **smaller number of interpretable patterns**, capturing regulatory programs or cell type–specific accessibility.  
- Mallet is a high-performance Java-based LDA implementation used here to efficiently handle large single-cell ATAC datasets.  


### 🧬 Mallet Example for Biologists

![Figure: MALLET NLP](MALLETNLP.png)
Source : Blei, D.M., 2012. Probabilistic topic models. Commun. ACM 55(4)


### In GRN context 



### 🧩 What is a GRN Module?

In the context of gene regulatory networks (GRNs), a **module** is:

> A group of **genes and/or regulatory regions (peaks)** that are **co-regulated**, meaning they tend to be active together in the same cells.  
> Modules often correspond to the set of genes controlled by the same transcription factor or regulatory program.  

Modules help simplify complex GRNs by grouping together genes with coordinated activity, making it easier to understand regulatory patterns in single-cell data.

| NLP Concept      | GRN Concept (SCENIC+) | What Mallet counts            | What LDA classifies             |
| ---------------- | --------------------- | ----------------------------- | ------------------------------- |
| Document         | Cell                  | Word occurrences per document | Topics in document              |
| Word             | Gene / Peak           | Frequency in document / cell  | Module membership               |
| Topic            | Module                | Co-occurring words/features   | Which genes/peaks form a module |
| Topic proportion | Module activity       | Counts in documents/cells     | Strength of module per cell     |



## Inputs

1. **Annotated cistopic object (`merged_with_meta.pkl`)**
   - Contains:
     - Peak-by-cell matrix (from previous steps)
     - scRNA-seq metadata for each cell (cell type, cluster, sample)

2. **Mallet software path**
   - The executable used to run LDA.

3. **Topic modeling parameters**
   - `n_topics`: number of topics to infer (e.g., 15, 20, 25, 30)
   - `n_iter`: number of iterations for the LDA algorithm
   - `alpha` / `eta`: Dirichlet priors controlling sparsity and topic distribution
   - `random_state`: ensures reproducibility

4. **Computational settings**
   - `n_cpu`: number of CPUs for parallel processing
   - `mallet_memory`: memory allocated for Mallet
   - Temporary and save directories (`tmp_path`, `save_path`)

---


## What it does?

1. **Prepare the peak-by-cell matrix**  
   - Converts the cistopic object data into a format compatible with Mallet.  

2. **Run LDA topic modeling**  
   - Mallet iteratively assigns peaks to topics based on co-accessibility patterns across cells.  
   - Each cell gets a **topic distribution** (proportion of reads assigned to each topic).  
   - Each topic gets a **peak distribution** (peaks enriched in that topic).  

3. **Store results**  
   - Saves Mallet output files in the specified directory, including topic assignments for cells and peaks.  

---

## Outputs

1. **Cistopic object with topic models**  
   - Contains updated metadata:
     - Each cell’s topic proportions  
     - Each topic’s peak composition  

2. **Mallet model files**  
   - Intermediate files for each topic and iteration  
   - Useful for diagnostics or rerunning analyses  

3. **Saved directories** (`MALLET/`)  
   - Organizes all output for downstream steps, including LDA visualization and DAR analysis.  

---

## Connection to Previous Steps

- **Takes the annotated cistopic object** from `add_scrna_metadata.py`  
- Uses the **peak-by-cell matrix** (from consensus peaks) and **cell metadata** (from scRNA annotations)  
- This is the core analytical step where the **high-dimensional ATAC-seq data is reduced into interpretable topics**, which will be used in clustering, DAR analysis, and visualization.

---

✅ **Summary:**  
- **Purpose:** Identify patterns of co-accessible genomic regions (topics) across cells.  
- **Inputs:** Annotated cistopic object, Mallet path, topic modeling parameters.  
- **Outputs:** Topic-assigned cistopic object, Mallet model files, ready for downstream analyses such as clustering and DAR detection.



# 11. Adding LDA Model to Cistopic Object

## Overview

This step integrates the **topic modeling results generated by Mallet** back into the cistopic object.  

- After running Mallet, the topic assignments and distributions are stored in external files.  
- To continue downstream pycisTopic analyses (clustering, DAR, visualization), the **cistopic object must include the LDA model results** internally.  


# 🧬 LDA in SCENIC+ (after Mallet)

**Input from Mallet:** A count matrix of genes/peaks per cell: each row = a cell, each column = a gene/peak, each entry = how often that gene/peak occurs (or is active) in that cell.

## Step 1: Initialize Modules (Topics)
LDA assumes that each cell is a mixture of modules (like a mixture of topics in NLP). It randomly assigns each gene/peak in each cell to a tentative module.  
Example:  
Cell 1: [GeneA, GeneB, GeneC, Peak1, Peak2] → [Module 0, Module 1, Module 0, Module 0, Module 1]  
Cell 2: [GeneD, GeneE, GeneC, Peak2, Peak3] → [Module 1, Module 1, Module 0, Module 1, Module 0]  
At this point, module assignments are just random guesses.

## Step 2: Iterative Update (Gibbs Sampling)
LDA repeatedly updates the assignment of each gene/peak to a module based on: 1) which modules are common in this cell (cell-level distribution) and 2) which modules usually contain this gene/peak across all cells (feature-level distribution). The algorithm uses probabilistic rules to favor assignments that explain the observed co-occurrences better. Intuition: Genes/peaks that tend to appear together across multiple cells gradually get grouped into the same module.

## Step 3: Define Modules
After many iterations, LDA converges. Each module is defined by a set of genes/peaks that co-occur across cells. These are the regulatory modules: groups of co-regulated genes/peaks.  
Example:  
Module 0: GeneA, GeneC, Peak1  
Module 1: GeneB, GeneD, GeneE, Peak2, Peak3

## Step 4: Module Activity per Cell
LDA calculates, for each cell, how strongly each module is represented:  
| Cell  | Module 0 | Module 1 |  
|-------|----------|----------|  
| 1     | 0.6      | 0.4      |  
| 2     | 0.4      | 0.6      |  
| 3     | 0.7      | 0.3      |  
This tells us which modules are active in each cell, analogous to topic proportions in NLP.

## ✅ Summary
Goal of LDA in SCENIC+: Discover hidden regulatory structure in single-cell data.  
Input: Mallet counts of genes/peaks per cell.  
Output: 1) Modules = groups of co-occurring/co-regulated genes/peaks, 2) Module activity per cell = how strongly each module is present in each cell. LDA is essentially finding hidden patterns of co-regulation from the raw co-occurrence data prepared by Mallet.



---

## Inputs

1. **Cistopic object with Mallet outputs (`MALLET/merged_cistopic_with_models.pkl`)**  
   - Contains:
     - Peak-by-cell matrix  
     - Cell metadata (including scRNA annotations)  
     - Placeholder for LDA topics  

2. **Output path (`outs/`)**  
   - Where the updated cistopic object with the integrated LDA model will be saved.

---

## What it does?

1. **Load the cistopic object and Mallet outputs**  
   - Reads the topic distributions for cells and peak compositions per topic.  

2. **Attach LDA results to the cistopic object**  
   - For each cell: stores **topic proportions**.  
   - For each topic: stores **peak weights / peak-to-topic assignment**.  

3. **Save the updated cistopic object**  
   - The object is now fully ready for:
     - Clustering based on topics  
     - Differential accessibility analysis (DAR)  
     - Visualization (UMAP, heatmaps, dotplots)  

---

## Outputs

1. **Cistopic object with integrated LDA model (`cistopic_LDA_obj.pkl`)**  
   - Each cell has topic proportions.  
   - Each topic has peak weights.  
   - Contains all previous metadata (scRNA annotations, QC info).  

---

## Connection to Previous Steps

- **Takes the cistopic object with merged cells and scRNA metadata**.  
- **Uses the topic modeling results generated by Mallet** in the previous step.  
- Prepares the object for **all downstream analyses in pycisTopic**, making it the “complete” cistopic object.

---

✅ **Summary:**  
- **Purpose:** Incorporate Mallet LDA results into the cistopic object for analysis.  
- **Inputs:** Cistopic object + Mallet topic results.  
- **Outputs:** Updated cistopic object with topics ready for clustering, DAR, and visualization.


# 12. Clustering Cistopic Objects and UMAP Visualization

## Overview

This step performs **dimensionality reduction and clustering** of cells based on their **topic profiles**.  

- Each cell now has a **topic proportion vector** from the previous LDA step.  
- Clustering identifies groups of cells with similar chromatin accessibility patterns.  
- UMAP is used to **visualize the cells in 2D space** for exploratory analysis and interpretation.  

---

## Inputs

1. **Cistopic object (`cistopic_obj_cleaned2.pkl`)**  
   - Contains:
     - Peak-by-cell matrix  
     - Cell metadata (including scRNA-seq annotations)  
     - LDA topic distributions per cell  

2. **Clustering parameters**  
   - `resolutions`: controls granularity of Leiden clustering (e.g., 0.6, 1.2, 3)  
   - `k`: number of neighbors used for graph construction  

3. **Output directories**  
   - Where UMAP plots and the clustered object will be saved.

---

## What it does?

1. **Dimensionality reduction**  
   - Computes a **cell-cell similarity graph** using the topic proportions.  

2. **Leiden clustering**  
   - Assigns cells into discrete clusters based on graph topology.  
   - Multiple resolutions allow exploring different cluster granularities.  

3. **UMAP embedding**  
   - Projects high-dimensional topic profiles into 2D for visualization.  
   - Each cell is positioned in 2D such that similar cells are close together.  

4. **Save results**  
   - Updates the cistopic object with **cluster assignments**.  
   - Saves UMAP coordinates and plots for exploration.

---

## Outputs

1. **Clustered cistopic object (`cistopic_obj_clustered.pkl`)**  
   - Each cell has:
     - Topic proportions  
     - Cluster assignment(s) at different resolutions  
     - Original metadata (scRNA annotations, QC info)  

2. **UMAP plots**  
   - Visualizes cell clusters and optionally sample or cell type labels.  

3. **UMAP coordinates and graph data**  
   - Stored in the object for downstream analyses.

---

## Connection to Previous Steps

- **Takes the cistopic object with integrated LDA model** (topic distributions per cell).  
- **Uses topic proportions** to construct a similarity graph for clustering.  
- Prepares the object for **differential accessibility analysis, visualization, and interpretation**.

---

✅ **Summary:**  
- **Purpose:** Group cells with similar chromatin accessibility patterns and visualize them in 2D.  
- **Inputs:** Cistopic object with topic proportions, clustering parameters.  
- **Outputs:** Clustered cistopic object with UMAP coordinates and plots, ready for downstream DAR analysis and visualization.

### Clustering output 

#### Cell level annotations
 
![metadata](outs/umap_clusters/metadata_umap.png?v=2)

![qc metrics](outs/umap_clusters/qc_metrics_umap.png?v=2)

#### Cluster level annotations 

![annotated clusters](outs/umap_clusters/annotated_clusters_umap.png?v=2)

![topics umap](outs/umap_clusters/topic_umap.png?v=2)

![topics heatmap ](outs/umap_clusters/cell_topic_heatmap.png?v=2)

#### Quick look to cell ratios in pycistopic 

![0.6](outs/stats/pycisTopic_leiden_10_0.6_composition_barplot.png)
![1.2](outs/stats/pycisTopic_leiden_10_1.2_composition_barplot.png)
![3.0](outs/stats/pycisTopic_leiden_10_3.0_composition_barplot.png)


# 13. Binarizing Topics Step

## Overview

This step **converts continuous topic distributions into binary accessibility matrices** for downstream analysis, such as differential accessibility testing (DAR).  

- Each cell has **topic proportions** from the LDA step.  
- Binarization transforms these proportions into **presence/absence calls**, e.g., a peak is considered "active" in a cell if its topic proportion passes a threshold.  
- This simplifies downstream analyses and makes them more robust to noise.

---

## Inputs

1. **Clustered cistopic object (`cistopic_obj_clustered.pkl`)**  
   - Contains:
     - Peak-by-cell matrix  
     - LDA topic proportions per cell  
     - Cluster assignments from the previous step  

2. **Output directory (`outs/topics`)**  
   - Where the binarized cistopic object and related files will be saved.

---

## What it does?

1. **Load clustered cistopic object**  
   - Access topic proportions for each cell.  

2. **Binarize topics**  
   - Convert continuous topic weights into binary values (0/1) for each peak in each cell.  
   - Thresholding ensures that only **high-confidence topic-peak associations** are considered active.  

3. **Update cistopic object**  
   - Replaces or adds a **binarized peak-by-cell matrix** to the object.  

4. **Save results**  
   - The binarized cistopic object is ready for **differential accessibility (DAR) analysis** and downstream analyses.

---

## Outputs

1. **Binarized cistopic object (`cistopic_obj_binarized.pkl`)**  
   - Contains binary peak-by-cell matrix alongside all previous metadata (topic assignments, clusters, scRNA annotations).  

2. **Optional summary files**  
   - Can include per-topic peak statistics or binarization thresholds.

---

## Connection to Previous Steps

- **Takes the clustered cistopic object** with topic proportions from the previous clustering step.  
- Converts continuous LDA results into a **binary format** suitable for DAR analysis.  
- Prepares the object for **robust differential accessibility testing** and other downstream applications.

---

✅ **Summary:**  
- **Purpose:** Convert continuous topic weights into binary activity calls for peaks.  
- **Inputs:** Clustered cistopic object with topic proportions.  
- **Outputs:** Binarized cistopic object ready for DAR analysis and downstream analyses.

#### Binarisations output 

![li](outs/topics/cell_topic_li.png?v=2)
![otsu](outs/topics/region_bin_otsu.png?v=2)
![top3k](outs/topics/region_bin_top3k.png?v=2)





# 14.  Differential Accessibility (DAR) Analysis Step

## Overview

This step identifies **differentially accessible regions (DARs)** between groups of cells, such as cell types or clusters.  

- Uses the **binarized cistopic object** from the previous step.  
- Tests which peaks are significantly more accessible in one group of cells versus others.  
- Results can reveal **cell type–specific regulatory regions** or other biologically meaningful patterns.  

---

## Inputs

1. **Binarized cistopic object (`cistopic_obj_binarized.pkl`)**  
   - Contains:
     - Binary peak-by-cell matrix  
     - Cell metadata (clusters, scRNA annotations, topic assignments)  

2. **Grouping variable (`-v celltype_scrna`)**  
   - Specifies which cell metadata column to use for group comparisons (e.g., cell type labels from scRNA-seq).  

3. **Analysis parameters**  
   - `n_cpu`: number of CPUs for parallel computation  
   - `temp_dir`: temporary files storage  
   - `scale_impute` / `scale_norm`: scaling factors for normalization  
   - `adjpval_thr`: adjusted p-value threshold for significance  
   - `log2fc_thr`: log2 fold-change threshold for significance  

4. **Output directory (`outs/DAR_results`)**  
   - Where DAR results and updated cistopic objects will be saved.

---

## What it does? 

1. **Load the binarized cistopic object**  
   - Access the binary peak-by-cell matrix and cell metadata.  

2. **Group cells by the specified variable**  
   - Example: group cells by `celltype_scrna`.  

3. **Compute differential accessibility**  
   - For each peak, compare accessibility across groups.  
   - Use statistical tests and thresholds (adj. p-value, log2 fold-change) to identify significant DARs.  

4. **Store results**  
   - Annotates the cistopic object with DAR information.  
   - Saves tables of significant DARs for each comparison.

---

## Outputs

1. **Cistopic object with DAR annotations (`cistopic_obj_with_DARs.pkl`)**  
   - Each peak is annotated with significance values, fold changes, and group-specific activity.  

2. **DAR result tables**  
   - Lists of peaks that are differentially accessible per group.  

3. **Optional plots or summaries**  
   - Can visualize DARs by cluster, cell type, or topic.

---

## Connection to Previous Steps

- **Takes the binarized cistopic object** from the previous binarization step.  
- Uses **cell type or cluster metadata** (from scRNA integration or clustering).  
- Generates **peak-level annotations** identifying regions differentially accessible between groups.  

---

✅ **Summary:**  
- **Purpose:** Detect peaks with significantly different accessibility across groups of cells.  
- **Inputs:** Binarized cistopic object, grouping variable, thresholds.  
- **Outputs:** Annotated cistopic object with DARs, result tables for downstream interpretation and visualization.


#### DAR preliminary results

![Higly variable Regions](outs/DAR_results/highly_variable_regions.png?v=4)

![Imputed features](outs/DAR_results/imputed_features.png?v=4)


# 15. Exporting Region Sets from DAR Results

## Overview

This step exports **lists of genomic regions (peaks) identified as DARs** into separate files for downstream analyses or external tools.  

- Takes the **cistopic object annotated with DAR results** from the previous step.  
- Generates **BED or other standard formats** representing differentially accessible regions per group or cell type.  
- These exported region sets can be used for:
  - Motif enrichment analysis  
  - Gene set enrichment analysis  
  - Visualization in genome browsers  

---

## Inputs

1. **Cistopic object with DAR annotations (`cistopic_obj_with_DARs.pkl`)**  
   - Contains:
     - Binary peak-by-cell matrix  
     - Cell metadata  
     - DAR results for each peak (p-values, fold-changes, significant groups)  

2. **Output directory (`outs/`)**  
   - Where the exported region sets will be saved.

---

## What it does?

1. **Load the DAR-annotated cistopic object**  
   - Access peaks and their associated DAR information.  

2. **Select peaks per group or significance criteria**  
   - Example: peaks with adjusted p-value below threshold and log2 fold-change above threshold.  

3. **Export to BED or other standard formats**  
   - One file per group or condition.  
   - Includes genomic coordinates and optionally additional metadata.  

4. **Save files in output directory**  
   - Ready for downstream analyses or sharing with other tools.

---

## Outputs

1. **Region set files (BED)**  
   - Lists of DARs per group, including genomic coordinates.  
   - Can be used for motif enrichment, visualization, or pathway analyses.  

2. **Optional summary tables**  
   - Summarize number of DARs per group or comparison.

---

## Connection to Previous Steps

- **Takes the DAR-annotated cistopic object** from `dar_analysis.py`.  
- Uses **peak-level differential accessibility information** to generate region sets.  
- Prepares results for **external analyses or downstream biological interpretation**.

---

✅ **Summary:**  
- **Purpose:** Export differentially accessible peaks as region sets for downstream analysis.  
- **Inputs:** Cistopic object with DAR annotations.  
- **Outputs:** BED or similar files containing DARs per group, ready for motif or pathway analyses.



# 🔴🔴🔴 PART C:  cisTarget Databases for Motif Enrichment (Optional Step)

## Overview

This step involves **using cisTarget databases** for motif and regulatory network analysis downstream of pycisTopic.  

- cisTarget databases link **genomic regions (peaks) to transcription factor motifs** and candidate target genes.  
- They are required for **motif enrichment analysis, regulatory network inference, and linking DARs to TF activity**.  
- **Creating custom cisTarget databases is optional**, but can be done if you want to analyze a specific genome, species, or motif collection.  
- For speed and convenience, we will **rely on prebuilt databases** in this workflow.

---

## Prebuilt Databases

For mouse (`mm10`) or human (`hg38`), prebuilt cisTarget databases typically include:

1. **Motif rankings (`.feather` files)**  
   - Genome-wide ranking of regions for each transcription factor motif.  
   - Example: `"mm10-500bp-upstream-7species.mc9nr.feather"`  

2. **Motif annotations (`.motifs.tbl`)**  
   - Information about the motifs, including TF name, family, and source database.  

3. **Gene annotations**  
   - Links between regions and nearby genes, allowing assignment of TF motifs to target genes.  

- These prebuilt files allow **rapid motif enrichment analysis** without having to scan the entire genome manually.

---

## Role in Workflow

- cisTarget databases are used **after DARs or topic regions are identified**.  
- They allow pycisTopic to perform:

  1. **Motif enrichment analysis**  
     - Determine which TF motifs are overrepresented in DARs or topic-specific peaks.  
  2. **Regulatory network inference**  
     - Identify candidate TFs controlling cell type–specific chromatin accessibility.  

- Without these databases, the downstream steps (cisTarget enrichment) cannot run efficiently.  

---

## Using Prebuilt vs Custom Databases

- **Prebuilt databases (recommended):**  
  - Saves time, standardized, widely used, compatible with pycisTopic.  
- **Custom databases (optional):**  
  - Can be created if you need specific motifs, a custom genome build, or updated annotations.  
  - Takes longer but provides maximum flexibility.

---

✅ **Summary:**  
- **Purpose:** Provide genome-wide motif rankings and annotations for enrichment analysis.  
- **Prebuilt inputs:** motif ranking files, motif annotations, and gene-region links.  
- **Outputs/usage:** Enable motif enrichment and TF-target network inference in downstream steps.  
- **Optional:** Users can build their own database for customized analyses.



🚀 **Note:** ===> I skipped this part and relied on prebuilt databases for speed 🚀




# 🔴🔴🔴 PART D: Running Scenic+ workflow step


## Overview

- The Snakemake workflow runs the SCENIC+ pipeline on your processed ATAC (and optional scRNA) data.  
- Its main goal is to **infer gene regulatory networks (GRNs)** and compute **TF activity scores per cell**.

---

## What it does?

1. **Prepare Input Matrices**
   - Extract peak-by-cell matrices or topic-specific binarized data from the cistopic object.  
   - Include cell metadata (clusters, cell types, scRNA info if available).  

2. **Motif Enrichment Using cisTarget**
   - For each DAR or topic-specific peak, identify enriched transcription factor motifs using prebuilt or custom cisTarget databases.  

3. **Regulatory Network Inference**
   - Combine motif enrichment with co-accessibility patterns (and optional scRNA co-expression) to link TFs to their predicted target genes.  
   - Construct a **gene regulatory network (GRN)**, identifying TF-gene relationships.  

4. **TF Activity Scoring**
   - Compute TF activity per cell based on accessibility of target peaks in the GRN.  
   - Produces a cell-by-TF activity matrix that can be used for clustering, visualization, or downstream analysis.  

5. **Visualization and Module Generation**
   - Generate regulatory modules: groups of co-regulated genes per TF.  
   - Create heatmaps, network diagrams, and motif enrichment plots for interpretation.

---

## Main Outputs

1. **Gene Regulatory Network (GRN)**
   - TFs linked to predicted target genes.  
   - Derived from DARs/topics and motif enrichment.  

2. **TF Activity Matrices**
   - Quantitative scores of TF activity per cell.  

3. **Regulatory Modules**
   - Groups of genes predicted to be co-regulated by each TF.  

4. **Plots and Visualizations**
   - Heatmaps of TF activity  
   - Network diagrams of GRNs  
   - Motif enrichment summaries  

---

## Connection to Previous Steps

- **Input:** binarized cistopic object with DARs or topic peaks (from pycisTopic).  
- **Input:** cisTarget databases (prebuilt or custom).  
- **Optional input:** scRNA-seq metadata for enhanced GRN inference.  
- **Output:** GRN, TF activity matrices, regulatory modules, and visualizations for downstream analysis.

---

✅ **Summary:**  
- **Purpose:** Infer gene regulatory networks and quantify TF activity per cell.  
- **Inputs:** Processed cistopic object (topics/DARs), cisTarget databases, optional scRNA metadata.  
- **Outputs:** GRN, TF activity matrices, regulatory modules, and visualizations.


## References

- [Pycistopic: Human Cerebellum Notebook](https://pycistopic.readthedocs.io/en/latest/notebooks/human_cerebellum.html#Getting-pseudobulk-profiles-from-cell-annotations)
- [SCENIC+: Official Documentation](https://scenicplus.readthedocs.io/en/latest/index.html)



