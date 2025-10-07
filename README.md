> # A walk through Scenic+ workflow: Gene Regulatory Network (GRN) Inference

### 🧠 Gene Regulatory Network (GRN) and Regulons

A **Gene Regulatory Network (GRN)** represents the regulatory relationships between **transcription factors (TFs)** and their **target genes**, describing how TFs can potentially control gene activity 

Within this network, a **regulon** is a functional unit consisting of a *TF and the set of genes it directly regulates*.  
In SCENIC+, regulons are identified and their activity is quantified across cells to reveal which regulatory programs are active in different cell states.

#### We will connect the pieces of the workflow, no optimization yet
  
## Installation hints 

```
snakemake --use-conda -j 1 -p \
    --resources mem_mb=256000 \
    --configfile Snakemake/config/config.yaml \
    --latency-wait 60

.snakemake/conda/
├── conda-abcdef123/   <- env for preprocessing.yml
├── conda-123456abc/   <- env for scenicplus.yml
├── conda-7890def456/   <- env for plotting.yml
├── conda-456def789/   <- env for some_other_rule.yml
├── conda-789abc012/   <- env for analysis_step2.yml
├── ...                <- env for visualization_extra.yml
├── ...
└── conda-zzz999xxx/   <- env for final_rule.yml

```

## 🔹 What SCENIC+ does

SCENIC+ is a workflow for inferring gene regulatory networks from single-cell multi-omics data, particularly scATAC-seq and multiome. It integrates chromatin accessibility with transcription factor motif analysis and gene expression information to identify *regulons* (transcription factors and their target genes). By combining topic modeling of accessibility profiles with motif enrichment and gene linkage, SCENIC+ maps regulatory programs to specific cell types and states, enabling the discovery of cell-type–specific transcriptional control.


> ### 🔹 Chromatin accessibility modeling (pycisTopic)
> - Perform topic modeling on scATAC-seq data to identify cis-regulatory topics (sets of co-accessible regions).
>
> ### 🔹 cisTarget database creation
> - Build motif/track databases from candidate regulatory regions.
> - Databases are used to connect TF motifs with accessible regions.
>
> ### 🔹 SCENIC+ analysis
> - Integrate ATAC-seq (regulatory regions) with RNA-seq (gene expression).
> - Infer TF–target gene regulons using motif enrichment + peak-to-gene links.
> - Score regulon activity across single cells to characterize cell states.

![Scenic+ workflow](scenicplusWorkflow.png)


---
# 🛑 PART A: Preprocessing of scRNA and scATAC
---

This step performs **cell-level preprocessing, dimensionality reduction, clustering, and visualization** for single-cell RNA-seq data.  

 
### 🔹 Role of scRNA-seq in SCENIC+

- **Linking regions to genes**  
  Uses gene expression to correlate accessible regions (from scATAC) with nearby genes, helping identify which open regions actually regulate those genes.
  And here we mean accessible regions are likely associated with gene activation not just close by the genome.

- **Scoring and validating regulatory modules**  
  Tests TF → region → gene modules against expression patterns — co-expressed target genes strengthen the TF–gene link.

- **Quantifying regulon activity per cell**  
  Computes regulon activity scores for each cell based on gene expression, enabling comparison across cell types or clusters to interpret biological function.


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


## scATAC preprocessing 

![ATAC UMAP](ATAC_samples.png)

### Clustering 

![ATAC CLUSTERS](ATAC_clusters.png)



---
# 🛑 PART B: Pycistopic
---


> ## 🔹 1. Pseudobulk Export with pycisTopic 

### 🔹 Initial Pseudobulk Aggregation (from fragments)

This step **combines raw ATAC-seq fragments into cell type–specific coverage tracks**:

- Reads from the **fragment files** of each sample are summed across all cells **.  
- Generates **BED/BigWig tracks** that capture chromatin accessibility at a summarized level.  

It uses `pycistopic_pseudobulk.py`

Total barcodes across all samples: 1156914
-----------------------
Sample distribution:
-----------------------
### 
| Sample   | Barcodes   |
|----------|------------|
| Control  | 584491     |
| KO       | 572423     |

and

```
scenicOuts/
└── consensus_peak_calling
    ├── bed_paths.tsv
    ├── bw_paths.tsv
    ├── pseudobulk_bed_files
    │   ├── Control.fragments.tsv.gz
    │   └── KO.fragments.tsv.gz
    └── pseudobulk_bw_files
        ├── Control.bw
        └── KO.bw
```
### A snapshot of the bigwig file 
![SCENIC+ BigWig Signal Tracks](chr3_filled_tracks.png)

> ## 🔹 2. Peak Calling Step with MACS2


This step identifies **peaks**, i.e., genomic regions that are significantly enriched for ATAC-seq fragments.  

- **Pseudobulk files** summarize chromatin accessibility for each **cell type × sample**.  
- **MACS2** runs **separately for each pseudobulk**, detecting regions with high accessibility that likely correspond to **regulatory elements** (enhancers or promoters).  
- The resulting **peak files (per sample)** are later combined to generate a **consensus peak set**, used for **topic modeling** in *pycisTopic*.


![Control](Control_peaks.png)  

![KO](KO_peaks.png)



> ## 🔹 3. Consensus Peak Generation Step

- Consensus peaks are created by merging all peaks detected across samples.  
- Any region open in **at least one sample** is included, forming a unified peak set.  
- This is useful when you want a shared reference space (same coordinates for all samples) for downstream analysis like LDA or accessibility matrices.
- Later analyses determine which peaks are condition-specific or shared.

---


> ## 🔹 4. TSS Generation Step in pycisTopic

This step generates a **BED file containing the transcription start sites (TSSs)** of genes for the reference genome.  

- TSS regions are important for **quality control**, such as checking **TSS enrichment** in ATAC-seq data.  
- It can also be used to **annotate peaks** with nearby genes for downstream analyses.  

> ## 🔹 5. QC check and plots 


### TH1
- Barcode QC  
  ![](scenicOuts/QC/TH1_barcode_qc.png)

- General QC  
  ![](scenicOuts/QC/TH1_qc.png)

### TH2
- Barcode QC  
  ![](scenicOuts/QC/TH2_barcode_qc.png)

- General QC  
  ![](scenicOuts/QC/TH2_qc.png)

## Comments on Graphs

| Plot | What to Look For | Good | Bad/Warning |
|------|-----------------|------|-------------|
| **TSS enrichment vs. unique fragments** | Relationship between signal strength and barcode count | Top-right cluster: strong enrichment, many fragments | Bottom-left: few fragments + low enrichment. Cluster at low enrichment but high fragments = possible junk |
| **FRiP vs. unique fragments** | Fraction of reads in peaks | High FRiP + enough fragments = high quality | Low FRiP or too few fragments = noisy |
| **Duplication ratio vs. unique fragments** | PCR duplication | Moderate duplication with enough fragments = ok | High duplication (PCR artifacts) or odd low-frag barcodes = junk |
| **Barcode rank plot** | Knee shape | Plateau (left side) = real cells | Steep drop (right) = background |
| **Fragment size distribution** | ATAC periodicity | Peaks at <100bp and ~200/400bp (nucleosomes) = good | Flat/no periodicity = poor complexity |
| **TSS profile** | Enrichment around TSS | Strong peak centered at 0 = good | Flat = background noise |


- A high number of fragments with a low fraction in peaks suggests that fragments are spread across large regions but aren't associated with meaningful or enriched features, indicating less informative data.
- High duplication ratios indicate that many fragments are duplicates, likely due to technical biases like PCR duplication, making the data less informative as they don't provide unique information.
- Before the Knee: The curve shows a sharp decline, indicating high-quality, informative barcodes with more unique and meaningful data.
- After the Knee: The curve flattens or drops sharply, suggesting low-quality data or background noise (e.g., technical artifacts, duplicates, low-coverage barcodes).


> ## 🔹 6. Creating Cistopic Objects Step

This step creates a **cistopic object**, which is the central data structure used by pycisTopic for **topic modeling of chromatin accessibility**.  

- The cistopic object organizes **fragment data, peak regions, and QC information** in a way suitable for downstream analyses.  
- It is essentially a **single-cell peak-by-cell matrix** stored in a Python pickle file, with metadata attached.  
- This step is **critical** because it transforms raw and pseudobulk fragment data into a structured object suitable for all downstream pycisTopic analyses.

## 🧫 Sample Overview

| Sample | Cells | Regions | Project Name | Fragment File |
|--------|-------|---------|--------------|---------------|
| **TH1** | 294 cells | 163,474 regions | `TH1` | `TH1_atac_fragments.tsv.gz` |
| **TH2** | 258 cells | 163,375 regions | `TH2` | `TH2_atac_fragments.tsv.gz` |


> ## 🔹 7. Merging Cistopic Objects Step

This step merges **one or more cistopic objects** into a single unified object.  

- In workflows with multiple samples, batches, or preprocessing runs, each cistopic object may represent a separate sample or subset of cells.  
- Merging combines them into a **single cistopic object**, making downstream analyses (topic modeling, clustering, DAR analysis) easier and consistent across all cells.  

## Merged CistopicObject Summary

| Metric          | Value                     |
|-----------------|---------------------------|
| **Cell Data**  | All 552 cells with complete QC metrics (23 columns) |
| **Total Regions** | 164,065                 |
| **Region Data**| Unified genomic regions from both samples       |


> ## 🔹 8. Adding scRNA-seq Metadata to Cistopic Objects

This step integrates **scRNA-seq-derived metadata** into the merged cistopic object.  

- Single-cell RNA-seq preprocessing (clustering, cell type annotation) provides **cell type labels, sample IDs, or other metadata**.  
- Attaching this information to the cistopic object allows **linking chromatin accessibility topics to known cell types** for interpretation.  


> ## 🔹 9. Topic Modeling with Mallet (run_mallet.py)

### 🧬 Latent Dirichlet Allocation (LDA) topic modeling/Mallet Example for Biologists

![Figure: LDA NLP](LDA.png)
Source : Blei, D.M., 2012. Probabilistic topic models. Commun. ACM 55(4)

Think of **documents** (e.g., research articles) as mixtures of **topics** (e.g., "neurodevelopment", "immune response").  
Each **topic** is defined by a group of **words that tend to occur together** (e.g., *axon, neuron, synapse*).

> ### 🔹 Latent Dirichlet Allocation (LDA) Topic Modeling
>
> 1. **🕵️ Guessing topics**  
>    - LDA scans across documents to find patterns of words that frequently co-occur.
>
> 2. **📝 Assigning words to topics**  
>    - Each word in a document is probabilistically assigned to a topic.  
>    - Assignments are refined using *Gibbs sampling* (shuffling topic labels until a stable pattern emerges).
>
> 3. **📊 Producing outputs**  
>    - A list of topics: groups of words that often appear together.  
>    - For each document: the proportion of topics it contains.

👉 In short: **LDA discovers hidden themes in large text collections by grouping words that “like to appear together,” and shows how much each document is made of those themes.**

---

### In GRN context 

```
> 🧩 `What is a GRN Module`
> A group of **genes and/or regulatory regions (peaks)** that are **co-regulated**, meaning they tend to be active together in the same cells.  
> Modules often correspond to the set of genes controlled by the same transcription factor or regulatory program.  
```

Modules help simplify complex GRNs by grouping together genes with coordinated activity, making it easier to understand regulatory patterns in single-cell data.

| NLP Concept      | GRN Concept (SCENIC+)       | What is counted (input to model)       | What LDA/MALLET infers                          |
| ---------------- | --------------------------- | -------------------------------------- | ------------------------------------------------ |
| Document         | Cell                        | Accessibility counts (peaks per cell) or expression counts (genes per cell) | Distribution of regulatory programs (topics) across the cell |
| Word             | Genomic region (peak) / Gene| How often a region is accessible in a cell (or gene is expressed) | Assignment of region/gene to one or more topics (modules) |
| Topic            | Regulatory module / Program | Groups of peaks or genes that co-occur across many cells | Which features form a coherent regulatory program |
| Topic proportion | Module activity in a cell   | Number of regions/genes linked to each topic in a cell | How strongly each regulatory program is active in that cell |

> 💡 **Note:** [MALLET](http://mallet.cs.umass.edu/) is a software package that provides an **efficient implementation of LDA**, making it faster and more scalable on large datasets.


``` 
python run_mallet.py \
  --cistopic_obj_pickle scenicOuts/merged_with_meta.pkl \
  --mallet_path /nfs/turbo/umms-thahoang/sherine/tools/Mallet-202108/bin/mallet \
  --n_topics 15 20 25 30 \
  --n_cpu 12 \
  --n_iter 500 \
  --tmp_path scenicOuts/TMP \
  --save_path scenicOuts/MALLET \
  --mallet_memory 300G \
  --random_state 555 \
  --alpha 5 \
  --alpha_by_topic \
  --eta 0.1 \
  --eta_by_topic
``` 

--- 
> ## 🔹 10. Adding LDA Model to Cistopic Object

**Input from Mallet:** A count matrix of genes/peaks per cell, where each row = a cell, each column = a gene/peak, and each entry = occurrence/activity.

**Goal of LDA in SCENIC+:** Discover hidden regulatory patterns in single-cell data.
- **Output:** 
  - 1) **Modules** = co-regulated gene/peak groups.
  - 2) **Module activity** = presence of each module in cells.

LDA identifies co-regulation patterns from Mallet's co-occurrence data.

### 🔹 LDA Topic Modeling on Chromatin Accessibility

LDA finds **patterns of coordinated accessibility** across cells:
- **Topics** represent groups of genomic regions (peaks) that open/close together.
- Each “topic” = one module of co-accessible peaks.
- For each cell, LDA measures the activity of each module.
  i.e. how much the cell shows the pattern of that module  
### 🔹 Output

Mallet (via LDA) generates two key matrices:
- **Region-Topic matrix:** Shows which peaks belong to each topic.
- **Topic-Cell matrix:** Shows how active each module is in each cell.

> ## 🔹 11. Clustering Cistopic Objects and UMAP Visualization

This step performs **dimensionality reduction and clustering** of cells based on their **topic profiles**.  

- Each cell now has a **topic proportion vector** from the previous LDA step.  
- Clustering identifies groups of cells with similar chromatin accessibility patterns.  
- UMAP is used to **visualize the cells in 2D space** for exploratory analysis and interpretation.  
- **Takes the cistopic object with integrated LDA model** (topic distributions per cell).
- **Uses topic proportions** to construct a similarity graph for clustering.
- Prepares the object for **differential accessibility analysis, visualization, and interpretation**.

### Clustering output 

###### Cell-based level 

![](scenicOuts/umap_clusters/celltype_umap.png)
###### Cluster-based level: More resolution more clusters more cell types as it is based on counting 
 
![](scenicOuts/umap_clusters/annotated_clusters_umap.png)
![](scenicOuts/umap_clusters/qc_metrics_umap.png)
![](scenicOuts/umap_clusters/topic_celltype_heatmap.png)


> ## 🔹 12. Binarizing Topics Step

This step **converts continuous topic distributions into binary accessibility matrices** for downstream analysis, such as differential accessibility testing (DAR).  

- Each cell has **topic proportions** from the LDA step.  
- Binarization transforms these proportions into **presence/absence calls**, e.g., a peak is considered "active" in a cell if its topic proportion passes a threshold.  
- This simplifies downstream analyses and makes them more robust to noise.

#### Binarisations output 

![Cell Topic LI](scenicOuts/topics/cell_topic_li.png)  


#### What the plots show

- **Each subplot** = one topic (e.g., Topic 1, Topic 2, … Topic 20).  
- **X-axis** = standardized probability of that topic in a given cell (range: 0–1).  
  - Higher probability → stronger association of the cell with the topic.  
- **Y-axis** = number of cells with that probability.  
  - Each subplot is a histogram of topic activity across cells.  
- **Red dashed vertical line** = binarisation threshold.  
  - Cells **to the right** of the line are considered *selected* (active for that topic).  
  - Cells **to the left** are not.  
- **“Selected = N”** = how many cells passed the threshold for that topic.  

---

#### How to interpret the plots

- **Topics with sharp right tails** (e.g., Topic 6, 13):  
  - Most cells have low probabilities; only a subset passes the threshold.  
  - → These topics are active in **specific subsets of cells** (likely distinct cell states or regulatory programs).  

- **Topics with broad or bimodal distributions** (e.g., Topic 4, 10):  
  - Many cells span medium to high probabilities.  
  - Threshold selects a large fraction of cells.  
  - → These topics may represent **widespread regulatory programs** across multiple cell types.  

- **Topics with flatter/noisy distributions** (e.g., Topic 2, 9):  
  - Fewer cells clearly pass the threshold.  
  - → These may be **weaker or less biologically meaningful topics**, or background noise.  

---

> ## 🔹 13.  Differential Accessibility (DAR) Analysis 

In this step, we try to find regions of the genome that are more accessible (open) in one group of cells compared to another. This helps us understand what makes different cell types or clusters unique.

- **Input:** We use a **binarized cistopic object** from the previous step (data that tells us which regions of the genome are open or closed in each cell).
  
- **Analysis:** We compare the accessibility of **genomic peaks** (regions) across different cell groups. The goal is to identify which peaks are more open in one group compared to another.

- **Output:** The results show which **genomic regions** are specific to certain cell types or groups, helping us uncover important biological differences between cells.

**Why it's different:** Unlike the previous step, which identifies general co-accessibility patterns across all cells, this step specifically compares accessibility across **different cell groups** to find regions that are uniquely accessible in one group versus others.


#### DAR preliminary results

![Higly variable Regions](scenicOuts/DAR_results/highly_variable_regions.png?v=4)

### Interpretation of DAR Mean–Dispersion Plot

* **Accessible regions:** The plot shows a large number of features, indicating plenty of accessible regions in the dataset.
* **Mean accessibility:** Most features have low mean accessibility (common in ATAC-seq), with some highly accessible peaks likely corresponding to promoters.
* **Dispersion:** A strong subset of features (yellow) are highly variable, suggesting good biological variability.
* **QC check:** The presence of many variable peaks (yellow) means the dataset has rich information for downstream DAR and regulatory network analysis.

✅ Overall: The data contain a healthy number of accessible and variable regions, suitable for SCENIC+ analysis.


![Imputed features](scenicOuts/DAR_results/imputed_features.png?v=4)

## Number of DARs found using  --adjpval_thr 0.05 and  --log2fc_thr 0.5

| Cell type   | DARs   |
|-------------|--------|
| AC          | 37,753 |
| BC          | 33,504 |
| Cones       | 0      |
| MG          | 23,660 |
| MGPC        | 13,113 |
| Microglia   | 15,838 |
| Rod         | 27,449 |


> ## 🔹 14. Exporting Region Sets from DAR Results

This step exports **lists of genomic regions (peaks) identified as DARs** into separate files for downstream analyses or external tools.  

- Takes the **cistopic object annotated with DAR results** from the previous step.  
- Generates **BED or other standard formats** representing differentially accessible regions per group or cell type.  
- These exported region sets can be used for:
  - Motif enrichment analysis  
  - Gene set enrichment analysis  
  - Visualization in genome browsers  

---
# 🛑 PART C: cisTarget Databases
---

- **cisTarget databases** link genomic regions (peaks) to **transcription factor (TF) motifs** **.  
They are essential for **motif enrichment analysis**, **regulatory network inference**, and for **linking DARs or accessible regions to potential TF activity**.

- **Creating custom cisTarget databases is optional**, but can be done if you want to analyze a specific genome, species, or motif collection.  


>> 1. 🧠 **The `.feather` file** tells you **“where might a motif bind?”** — i.e., which genomic regions look like good binding sites for each motif (based purely on DNA sequence similarity).

>> 2. 🧠 **The `.motifs.tbl` file** tells you **“which transcription factor probably uses that motif?”** — i.e., it adds biological context by mapping motifs to their likely TF(s), similar motifs, and orthologs.


### Feather file Example

| chr10:100001588-100001754 | chr10:100004590-100004749 | chr10:100009729-100010013 | ... | chrY:9986285-9986625 | chrY:9992323-9992533 | motifs          |
|----------------------------|---------------------------|---------------------------|-----|----------------------|----------------------|-----------------|
| 635378                     | 961923                    | 81883                     | ... | 261962               | 1096188              | bergman__Su_H_  |
| 294685                     | 79408                     | 809215                    | ... | 108758               | 747006               | bergman__croc   |
| 577456                     | 435983                    | 23497                     | ... | 682309               | 706782               | bergman__tll    |
| 976643                     | 795153                    | 810566                    | ... | 3809                 | 1018969              | c2h2_zfs__M0369 |


### tbl file Example 
 
| motif_id         | motif_name | gene_name | motif_similarity_qvalue | similar_motif_id | orthologous_identity |
| ---------------- | ---------- | --------- | ----------------------- | ---------------- | -------------------- |
| metacluster_196.3| EcR_usp    | Hnf4a     | 0                       | None             | 0.265823             |
| metacluster_196.3| EcR_usp    | Nr1h4     | 1.04865e-07             | metacluster_64.14| 0.539514             |
| metacluster_196.3| EcR_usp    | Nr2f1     | 8.38436e-08             | metacluster_64.14| 0.995238             |


---
# 🛑 PART D: Running Scenic+ workflow step 
---

SCENIC+ builds on the outputs of **pycistopic** and **cistarget** to generate gene regulatory networks (GRNs).
----

> **Outputs**  
>  
> 1. **Gene Regulatory Network (GRN)**  
>    - TFs linked to predicted target genes  
>    - Derived from DARs/topics and motif enrichment  
>  
> 2. **TF Activity Matrices**  
>    - Quantitative scores of TF activity per cell  
>  
> 3. **Regulatory Modules**  
>    - Groups of genes predicted to be co-regulated by each TF  
>  
> 4. **Plots and Visualizations**  
>    - Heatmaps of TF activity  
>    - Network diagrams of GRNs  
>    - Motif enrichment summaries  

----

#### What SCENIC+ does

👉 In short, SCENIC+ moves from **chromatin accessibility → motif enrichment → region-to-gene linking → TF activity per cell**, providing a powerful way to dissect **cell type–specific gene regulation**.  


#### 1. Motif Enrichment (cisTarget step)

🔹 **cisTarget (ctx):**  
	Finds enriched motifs
        ![](ctx_sample.png)
🔹 **Differential Motif Enrichment (DEM):**  
       Finds motifs that change activity between groups.
       ![](dem_sample.png)
        

📄 The **ctx output** — [`ctx_results.html`](Snakemake/workflow/ctx_results.html) —  
summarizes **motif enrichment** across region sets (e.g., topics or DARs), helping identify transcription factors whose motifs are overrepresented in accessible regions.

📄 The **DEM output** — [`dem_results.html`](Snakemake/workflow/dem_results.html) —  
reports **differential motif activity** across cell groups or conditions, highlighting motifs with cluster- or condition-specific accessibility patterns.

---

#### 2. Link Regions to Genes
- Connect enriched regions to their **nearby or co-expressed genes**, using both **genomic proximity** and **correlation with RNA expression**.  
- Produces **region–gene relationships** that suggest potential regulatory links.

###### `region_to_gene` Example

| target        | region                 |           importance |                  rho |       importance_x_rho |  importance_x_abs_rho | Distance  |
| ------------- | ---------------------- | -------------------: | -------------------: | ---------------------: | --------------------: | --------- |
| 0610005C13Rik | chr7:45451023-45451660 | 0.023644760576675295 | -0.07113737037272991 | -0.0016820260905174732 | 0.0016820260905174732 | [-116452] |
| 0610005C13Rik | chr7:45588873-45589602 | 0.017281679828869717 |  0.08660415959528277 |  0.0014966653579740119 | 0.0014966653579740119 | [14062]   |
| 0610005C13Rik | chr7:45570041-45570645 | 0.016397470259162692 |  0.10171570471952794 |  0.0016678802430282337 | 0.0016678802430282337 | [0]       |
| 0610005C13Rik | chr7:45567744-45568545 | 0.019381896115026935 |  0.05951350104480214 |   0.001153484494691902 |  0.001153484494691902 | [0]       |

**Description:** Links ATAC peaks (regions) to the target gene **0610005C13Rik**.

* **importance:** model-derived weight (how informative the region is for predicting the gene).
* **rho:** Spearman correlation between region accessibility and gene expression (sign shows direction).
* **importance_x_rho:** signed score combining weight and correlation (positive → consistent with activation).
* **importance_x_abs_rho:** absolute-strength score (magnitude only).
* **Distance:** distance from region to the gene TSS (bp); negative = upstream, positive = downstream.
* “Distance = 0” = region overlaps the gene’s promoter (very likely a direct regulatory link). 
---

### ⏳ In Progress / Pending 

#### 3. Build Regulatory Networks
- Combine **TF–region links** (from motif enrichment) with **region–gene links**.  
- Result: TF → region → gene connections, i.e. a **Gene Regulatory Network (GRN)**.  
- Optionally integrate **scRNA-seq expression** to refine TF–target predictions.

---

#### 4. Score TF Activity Per Cell
- Evaluate the accessibility of each TF’s target regions per cell.  
- Produces a **cell × TF activity matrix** (similar to regulon activity in SCENIC).  
- Allows clustering and visualization of regulatory programs across cell states.

---

#### 5. Visualization & Modules
- Summarize TF–target relationships into **regulatory modules** (groups of co-regulated genes).  
- Generate plots:  
  - Heatmaps of TF activity  
  - Motif enrichment plots  
  - Network diagrams of TFs and targets  

---


## References

- [Pycistopic: Human Cerebellum Notebook](https://pycistopic.readthedocs.io/en/latest/notebooks/human_cerebellum.html#Getting-pseudobulk-profiles-from-cell-annotations)
- [SCENIC+: Official Documentation](https://scenicplus.readthedocs.io/en/latest/index.html)



