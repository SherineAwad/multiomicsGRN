> # A walk through Scenic+ workflow: Gene Regulatory Network (GRN) Inference

#### ⚠️ Disclaimer: We will connect the pieces of the workflow, no optimization yet

### 🧠 Gene Regulatory Network (GRN) and Regulons

A **Gene Regulatory Network (GRN)** represents the regulatory relationships between **transcription factors (TFs)** and their **target genes**, describing how TFs can potentially control gene activity 

Within this network, a **regulon** is a functional unit consisting of a *TF and the set of genes it directly regulates*.  
In SCENIC+, regulons are identified and their activity is quantified across cells to reveal which regulatory programs are active in different cell states.

#### We will connect the pieces of the workflow, no optimization yet
  
## Workflow -- Step A-Z conda envs  

```
snakemake --use-conda -j 1 -p \
    --resources mem_mb=256000 \
    --configfile Snakemake/config/config.yaml \
    --latency-wait 60

.snakemake/conda/
├── conda-abcdef123/   <-  env1.yml env for step1,2, and 3 
├── conda-123456abc/   <-  env2.yml env for step4,5, and 6 
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

##  scRNA annotations: to use a h5ad accepted by scenic+ 

![Clusters UMAP](figures/umap_clustered_sNeurog2_Clusters.png)

### DOTPLOT 
![Dotplot](figures/clustered_sNeurog2_Dotplot.png)

### Marker Genes UMAP

<img src="figures/umap_clustered_sNeurog2_Abca8a.png?v=4" alt="Abca8a" width="33%"><img src="figures/umap_clustered_sNeurog2_Ccr2.png?v=4" alt="Ccr2" width="33%"><img src="figures/umap_clustered_sNeurog2_Hes1.png?v=4" alt="Hes1" width="33%">
<img src="figures/umap_clustered_sNeurog2_Notch1.png?v=4" alt="Notch1" width="33%"><img src="figures/umap_clustered_sNeurog2_Rlbp1.png?v=4" alt="Rlbp1" width="33%"><img src="figures/umap_clustered_sNeurog2_Acta2.png?v=4" alt="Acta2" width="33%">
<img src="figures/umap_clustered_sNeurog2_Chat.png?v=4" alt="Chat" width="33%"><img src="figures/umap_clustered_sNeurog2_Hes5.png?v=4" alt="Hes5" width="33%"><img src="figures/umap_clustered_sNeurog2_Nrl.png?v=4" alt="Nrl" width="33%">
<img src="figures/umap_clustered_sNeurog2_Rpe65.png?v=4" alt="Rpe65" width="33%"><img src="figures/umap_clustered_sNeurog2_Apoe.png?v=4" alt="Apoe" width="33%">
<img src="figures/umap_clustered_sNeurog2_Insm1.png?v=4" alt="Insm1" width="33%"><img src="figures/umap_clustered_sNeurog2_Olig2.png?v=4" alt="Olig2" width="33%"><img src="figures/umap_clustered_sNeurog2_Sebox.png?v=4" alt="Sebox" width="33%">
<img src="figures/umap_clustered_sNeurog2_Aqp4.png?v=4" alt="Aqp4" width="33%"><img src="figures/umap_clustered_sNeurog2_Csf1r.png?v=4" alt="Csf1r" width="33%"><img src="figures/umap_clustered_sNeurog2_Isl1.png?v=4" alt="Isl1" width="33%">
<img src="figures/umap_clustered_sNeurog2_Otx2.png?v=4" alt="Otx2" width="33%"><img src="figures/umap_clustered_sNeurog2_Slc17a7.png?v=4" alt="Slc17a7" width="33%"><img src="figures/umap_clustered_sNeurog2_Arr3.png?v=4" alt="Arr3" width="33%">
<img src="figures/umap_clustered_sNeurog2_Elavl3.png?v=4" alt="Elavl3" width="33%"><img src="figures/umap_clustered_sNeurog2_Kcnj8.png?v=4" alt="Kcnj8" width="33%"><img src="figures/umap_clustered_sNeurog2_Pax2.png?v=4" alt="Pax2" width="33%">
<img src="figures/umap_clustered_sNeurog2_Slc1a3.png?v=4" alt="Slc1a3" width="33%"><img src="figures/umap_clustered_sNeurog2_Ascl1.png?v=4" alt="Ascl1" width="33%"><img src="figures/umap_clustered_sNeurog2_Elavl4.png?v=4" alt="Elavl4" width="33%">
<img src="figures/umap_clustered_sNeurog2_Lhx1.png?v=4" alt="Lhx1" width="33%"><img src="figures/umap_clustered_sNeurog2_Pax6.png?v=4" alt="Pax6" width="33%"><img src="figures/umap_clustered_sNeurog2_Slc6a9.png?v=4" alt="Slc6a9" width="33%">
<img src="figures/umap_clustered_sNeurog2_Atoh7.png?v=4" alt="Atoh7" width="33%"><img src="figures/umap_clustered_sNeurog2_Emx1.png?v=4" alt="Emx1" width="33%"><img src="figures/umap_clustered_sNeurog2_Lhx2.png?v=4" alt="Lhx2" width="33%">
<img src="figures/umap_clustered_sNeurog2_Pou4f2.png?v=4" alt="Pou4f2" width="33%"><img src="figures/umap_clustered_sNeurog2_Sox11.png?v=4" alt="Sox11" width="33%"><img src="figures/umap_clustered_sNeurog2_Bsn.png?v=4" alt="Bsn" width="33%">
<img src="figures/umap_clustered_sNeurog2_Foxn4.png?v=4" alt="Foxn4" width="33%"><img src="figures/umap_clustered_sNeurog2_Lhx4.png?v=4" alt="Lhx4" width="33%"><img src="figures/umap_clustered_sNeurog2_Prdm1.png?v=4" alt="Prdm1" width="33%">
<img src="figures/umap_clustered_sNeurog2_Sox9.png?v=4" alt="Sox9" width="33%"><img src="figures/umap_clustered_sNeurog2_Cabp5.png?v=4" alt="Cabp5" width="33%"><img src="figures/umap_clustered_sNeurog2_Gad1.png?v=4" alt="Gad1" width="33%">
<img src="figures/umap_clustered_sNeurog2_Malat1.png?v=4" alt="Malat1" width="33%"><img src="figures/umap_clustered_sNeurog2_Prdx6.png?v=4" alt="Prdx6" width="33%"><img src="figures/umap_clustered_sNeurog2_Tfap2a.png?v=4" alt="Tfap2a" width="33%">
<img src="figures/umap_clustered_sNeurog2_Calb1.png?v=4" alt="Calb1" width="33%"><img src="figures/umap_clustered_sNeurog2_Gfap.png?v=4" alt="Gfap" width="33%"><img src="figures/umap_clustered_sNeurog2_mt-Atp6.png?v=4" alt="mt-Atp6" width="33%">
<img src="figures/umap_clustered_sNeurog2_Rbfox3.png?v=4" alt="Rbfox3" width="33%"><img src="figures/umap_clustered_sNeurog2_Tie1.png?v=4" alt="Tie1" width="33%"><img src="figures/umap_clustered_sNeurog2_Calb2.png?v=4" alt="Calb2" width="33%">
<img src="figures/umap_clustered_sNeurog2_Glul.png?v=4" alt="Glul" width="33%"><img src="figures/umap_clustered_sNeurog2_Neurog2.png?v=4" alt="Neurog2" width="33%"><img src="figures/umap_clustered_sNeurog2_Rho.png?v=4" alt="Rho" width="33%">
<img src="figures/umap_clustered_sNeurog2_Vim.png?v=4" alt="Vim" width="33%">

### Annotations 

![Annotations ON](figures/annotated_clustered_sNeurog2_annotationsON.png?v=3)

![Annotations](figures/annotated_clustered_sNeurog2_annotations.png?v=3)

##### 🧬 Cell Counts per Type a
| Cell Type   | Count |
|-------------|--------|
| MG          | 10043  |
| BC          | 3681   |
| MGPC        | 1662   |
| Rod         | 1385   |
| Microglia   | 1189   |
| AC          | 571    |
| Cones       | 95     |



##### 🧬 Cell Counts per Typpe and Sample 

| Cell Type   | Control Count | OX Count |
|-------------|---------------|----------|
| MGPC        | 8             | 1654     |
| AC          | 27            | 544      |
| BC          | 132           | 3549     |
| Cones       | 13            | 82       |
| MG          | 7419          | 2624     |
| Microglia   | 957           | 232      |
| Rod         | 505           | 880      |




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
scenicResults200/
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

##### Chr1 [--start 97,710,985 --end 97,760,985 ]
![chr1](chr1_middle.png)

##### Chr19 
![chr19](chr19_full.png)

> ## 🔹 2. Peak Calling Step with MACS2


This step identifies **peaks**, i.e., genomic regions that are significantly enriched for ATAC-seq fragments.  

- **Pseudobulk files** summarize chromatin accessibility for each **cell type × sample**.  
- **MACS2** runs **separately for each pseudobulk**, detecting regions with high accessibility that likely correspond to **regulatory elements** (enhancers or promoters).  
- The resulting **peak files (per sample)** are later combined to generate a **consensus peak set**, used for **topic modeling** in *pycisTopic*.

###### Control Peaks 
![Control](Control_peaks.png)  
###### KO Peaks
![KO](KO_peaks.png)



> ## 🔹 3. Consensus Peak Generation Step

- Peaks from all samples are merged to form a **unified peak set** representing all accessible regions.
- Any region open in **at least one sample** is included, ensuring no regulatory element is missed.
- Overlapping peaks are resolved by:
  - Keeping the original peak if only one overlaps,
  - Selecting the peak with the highest score if two overlap,
  - Choosing the most significant peak and removing overlapping others if three or more overlap.
- This creates a **shared reference space** (same peak coordinates across samples) for downstream analyses like LDA and accessibility matrices.
- Later steps identify which peaks are condition-specific or shared across samples.

---


> ## 🔹 4. TSS Generation Step in pycisTopic

This step generates a **BED file containing the transcription start sites (TSSs)** of genes for the reference genome.  

- TSS regions are important for **quality control**, such as checking **TSS enrichment** in ATAC-seq data.  
- It can also be used to **annotate peaks** with nearby genes for downstream analyses.  

> ## 🔹 5. QC check and plots 

```
python collect_qc_barcodes.py \
    --fragments_dict scenicResults200/fragments_dict.pkl \
    --qc_output_dir scenicResults200/QC \
    --output_pickle scenicResults200/QC/qc_barcodes_thresholds.pkl \
    --unique_fragments_threshold 200 \
    --tss_enrichment_threshold 0.1 \
    --frip_threshold 0
```

### TH1
- Barcode QC  
  ![](scenicResults200/QC/TH1_barcode_qc.png)

- General QC  
  ![](scenicResults200/QC/TH1_qc.png)

### TH2
- Barcode QC  
  ![](scenicResults200/QC/TH2_barcode_qc.png)

- General QC  
  ![](scenicResults200/QC/TH2_qc.png)

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

###### How PycisTopic Calculates the `unique_fragments_threshold`

- PycisTopic determines the `unique_fragments_threshold` **per sample** based on the distribution of fragments per barcode.
- The process roughly works as follows:
  1. **Count unique fragments per barcode:** For each barcode in the ATAC-seq fragment file, count the number of distinct fragments mapped to the genome.
  2. **Estimate a quality threshold:** Using the distribution of fragment counts, PycisTopic selects a threshold to filter out low-quality barcodes (likely empty droplets or background).
     - This can be percentile-based or use heuristics on the fragment count distribution.
  3. **Apply the threshold:** Barcodes with **unique fragments below this value** are removed from the dataset.



> ## 🔹 6. Creating Cistopic Objects Step

This step creates a **cistopic object**, which is the central data structure used by pycisTopic for **topic modeling of chromatin accessibility**.  

- The cistopic object organizes **fragment data, peak regions, and QC information** in a way suitable for downstream analyses.  
- It is essentially a **single-cell peak-by-cell matrix** stored in a Python pickle file, with metadata attached.  
- This step is **critical** because it transforms raw and pseudobulk fragment data into a structured object suitable for all downstream pycisTopic analyses.


> ## 🔹 7. Merging Cistopic Objects Step

This step merges **one or more cistopic objects** into a single unified object.  

- In workflows with multiple samples, batches, or preprocessing runs, each cistopic object may represent a separate sample or subset of cells.  
- Merging combines them into a **single cistopic object**, making downstream analyses (topic modeling, clustering, DAR analysis) easier and consistent across all cells.  


#### 📊 Sample Overview

- **Sample 1**: 7,517 cells  
- **Sample 2**: 7,675 cells  

#### 📋 Summary Statistics

- **Total number of cells**: 15,192  
- **Total number of regions**: 164,077



> ## 🔹 8. Adding scRNA-seq Metadata to Cistopic Objects

This step integrates **scRNA-seq-derived metadata** into the merged cistopic object.  

- Single-cell RNA-seq preprocessing (clustering, cell type annotation) provides **cell type labels, sample IDs, or other metadata**.  
- Attaching this information to the cistopic object allows **linking chromatin accessibility topics to known cell types** for interpretation.  

#### After merging 

#### 📊 Sample Overview

- **Sample 1**: 6,537 cells  
- **Sample 2**: 6,542 cells  

#### 📋 Summary Statistics

- **Total number of cells**: 13,079 (86.05% match) 
- **Total number of regions**: 164,077


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
> A group of **genes and/or regulatory regions (peaks)** that are **co-regulated**, 
  meaning they tend to be active together in the same cells.  
> Modules often correspond to the set of genes controlled by 
  the same transcription factor or regulatory program.  
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
  --cistopic_obj_pickle scenicResults200/merged_with_meta.pkl \
  --mallet_path /nfs/turbo/umms-thahoang/sherine/tools/Mallet-202108/bin/mallet \
  --n_topics 5 15 25 45 55 \
  --n_cpu 12 \
  --n_iter 500 \
  --tmp_path scenicResults200/TMP \
  --save_path scenicResults200/MALLET \
  --mallet_memory 300G \
  --random_state 555 \
  --alpha 1 \
  --alpha_by_topic \
  --eta 0.1 \
  --eta_by_topic
``` 

--- 
> ## 🔹 10. Adding LDA Model to Cistopic Object

LDA finds **patterns of coordinated accessibility** across cells:
- **Topics** represent groups of genomic regions (peaks) that open/close together.
- Each “topic” = one module of co-accessible peaks.
- For each cell, LDA measures the activity of each module.
  i.e. how much the cell shows the pattern of that module  

> ## 🔹 11. Clustering Cistopic Objects and UMAP Visualization

This step performs **dimensionality reduction and clustering** of cells based on their **topic profiles**.  

- Each cell now has a **topic proportion vector** from the previous LDA step.  
- Clustering identifies groups of cells with similar chromatin accessibZZility patterns.  

### Clustering output 

###### Cell-based level 

![](scenicResults200/umap_clusters/celltype_umap.png)
 
![](scenicResults200/umap_clusters/topic_celltype_heatmap.png) 
![](scenicResults200/umap_clusters/annotated_clusters_res_3.0.png)
![](scenicResults200/umap_clusters/annotated_clusters_res_2.0.png)
![](scenicResults200/umap_clusters/annotated_clusters_res_1.5.png)
![](scenicResults200/umap_clusters/annotated_clusters_res_1.0.png)
![](scenicResults200/umap_clusters/annotated_clusters_res_0.8.png)
![](scenicResults200/umap_clusters/annotated_clusters_res_0.5.png)


> ## 🔹 12. Binarizing Topics Step

This step **converts continuous topic distributions into binary accessibility matrices** for downstream analysis, such as differential accessibility testing (DAR).  

- Each cell has **topic proportions** from the LDA step.  
- Binarization transforms these proportions into **presence/absence calls**, e.g., a peak is considered "active" in a cell if its topic proportion passes a threshold.  
- This simplifies downstream analyses and makes them more robust to noise.

#### Binarisations output 

![Cell Topic LI](scenicResults200/topics/cell_topic_li.png?v=4)

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

- **Topics with sharp right tails** (e.g., Topic 16):  
  - Most cells have low probabilities; only a subset passes the threshold.  
  - → These topics are active in **specific subsets of cells** (likely distinct cell states or regulatory programs).  

- **Topics with broad or bimodal distributions** (e.g., Topic 2):  
  - Many cells span medium to high probabilities.  
  - Threshold selects a large fraction of cells.  
  - → These topics may represent **widespread regulatory programs** across multiple cell types.  

- **Topics with flatter/noisy distributions** (e.g., Topic 12):  
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

## Number of DARs found using  --adjpval_thr 0.05 and  --log2fc_thr 0.5

## Number of DARs found

| Cell Type  | Number of DARs |
|------------|----------------|
| AC         | 14,986         |
| BC         | 35,141         |
| Cones      | 28,307         |
| MG         | 19,906         |
| MGPC       | 8,535          |
| Microglia  | 18,403         |
| Rod        | 21,095         |


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
| 6 35378                     | 961923                    | 81883                     | ... | 261962               | 1096188              | bergman__Su_H_  |
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

### What SCENIC+ does


```plaintext
Chromatin Accessibility  -->  Motif Enrichment  -->  Region-to-Gene Linking  -->  TF Activity per Cell
```

##### Full Snakemake DAG 

![](snake_dag.png)


### 1. Motif Enrichment (cisTarget step)

🔹 **cisTarget (ctx):**  
        Find motifs enriched in accessible regions
        ![](ctx_sample.png)

   - **NES (Normalized Enrichment Score):** How strongly the motif is enriched compared to background — higher = stronger enrichment.
   - **AUC (Area Under Curve):** Reflects how consistently the motif ranks high across regions — measures enrichment quality.
   - **Rank_at_max:** The position in the ranked region list where enrichment was highest.

🔹 **Differential Motif Enrichment (DEM):**  
      Finds motifs with activity differences between groups or conditions. 
      ![](dem_sample.png)
        

📄 The **ctx output** 
summarizes **motif enrichment** across region sets (e.g., topics or DARs), helping identify transcription factors whose motifs are overrepresented in accessible regions.

📄 The **DEM output** 
reports **differential motif activity** across cell groups or conditions, highlighting motifs with cluster- or condition-specific accessibility patterns.


## 🐛 Debugging / Logs

### Example Regions
These are some of the genomic regions being analyzed:

| Region |
|--------|
| chr10:100009435-100010639 |
| chr10:100015388-100016745 |
| chr10:100017719-100018025 |
| chr10:100021901-100022856 |

### Feather File Format
The feather file contains **region counts** per genomic bin (integer values):

| Region | Data Type |
|--------|-----------|
| chr10:100001588-100001754 | int32 |
| chr10:100004590-100004749 | int32 |
| chr10:100009729-100010013 | int32 |
| chr10:100010349-100010586 | int32 |

**✅ Note:** This format is compatible with SCENIC+ and can be used for region-to-gene inference.

= 228 TFs detected

---

### 2. Link Regions to Genes
- Connect enriched regions to their **nearby or co-expressed genes**, using both **genomic proximity** and **correlation with RNA expression**.  
- Produces **region–gene relationships** that suggest potential regulatory links.
- “These peaks have binding motifs for TF X — which genes might they be regulating?” 
- Assign the region to nearby genes (e.g., within 10–100 kb) 
- Check if accessibility of the region is correlated with the expression of a gene across single cells

###### `region_to_gene` Example

| target        | region                 |           importance |                  rho |       importance_x_rho |  importance_x_abs_rho | Distance  |
| ------------- | ---------------------- | -------------------: | -------------------: | ---------------------: | --------------------: | --------- |
| 0610005C13Rik | chr7:45451023-45451660 | 0.023644760576675295 | -0.07113737037272991 | -0.0016820260905174732 | 0.0016820260905174732 | [-116452] |
| 0610005C13Rik | chr7:45588873-45589602 | 0.017281679828869717 |  0.08660415959528277 |  0.0014966653579740119 | 0.0014966653579740119 | [14062]   |
| 0610005C13Rik | chr7:45570041-45570645 | 0.016397470259162692 |  0.10171570471952794 |  0.0016678802430282337 | 0.0016678802430282337 | [0]       |
| 0610005C13Rik | chr7:45567744-45568545 | 0.019381896115026935 |  0.05951350104480214 |   0.001153484494691902 |  0.001153484494691902 | [0]       |

**Description:** Links ATAC peaks (regions) to the target gene **0610005C13Rik**.

* **importance:** model-derived weight - How important the motif/TF binding in this region is (from earlier motif enrichment step) 
* **rho:** Spearman correlation between region accessibility and gene expression (sign shows direction).
* **importance_x_rho:** signed score combining weight and correlation (positive → consistent with activation).
* **importance_x_abs_rho:** absolute-strength score (magnitude only).
* **Distance:** distance from region to the gene TSS (bp); negative = upstream, positive = downstream.
* “Distance = 0” = region overlaps the gene’s promoter (very likely a direct regulatory link). 
---

### 3. Build Regulatory Networks
- Combine **TF–region links** (from motif enrichment) with **region–gene links**.  
- Result: TF → region → gene connections, i.e. a **Gene Regulatory Network (GRN)**.  
- Optionally integrate **scRNA-seq expression** to refine TF–target predictions.

---

### 4. Score TF Activity Per Cell
- Evaluate the accessibility of each TF’s target regions per cell.  
- Produces a **cell × TF activity matrix** (similar to regulon activity in SCENIC).  
- Allows clustering and visualization of regulatory programs across cell states.

---

### 5. Visualization & Modules
- Summarize TF–target relationships into **regulatory modules** (groups of co-regulated genes).  
- Generate plots:  
  - Heatmaps of TF activity  
  - Motif enrichment plots  
  - Network diagrams of TFs and targets  

--- 

#### Snakefile rules

| Job                             | Count |
|---------------------------------|-------|
| AUCell_direct                    | 1     |
| AUCell_extended                  | 1     |
| all                              | 1     |
| eGRN_direct                      | 1     |
| eGRN_extended                    | 1     |
| get_search_space                 | 1     |
| motif_enrichment_cistarget       | 1     |
| motif_enrichment_dem             | 1     |
| prepare_GEX_ACC_non_multiome     | 1     |
| prepare_menr                     | 1     |
| region_to_gene                   | 1     |
| scplus_mudata                    | 1     |
| tf_to_gene                       | 1     |
| total                            | 13    |


## First output: eRegulon_Direct.tsv

### SCENIC+ `eRegulon_direct.tsv` 

The `eRegulon_direct.tsv` file captures **direct regulatory links** between genomic regions, transcription factors (TFs), and their target genes. Each row represents a **triplet**:  

**Region → TF → Gene**

---

## Column Descriptions

| Column | Meaning |
|--------|---------|
| **Region** | Genomic coordinates of the regulatory element (`chr:start-end`). |
| **Gene** | Target gene associated with the region. |
| **importance_R2G** | Importance of the **Region → Gene** link (higher = stronger influence). |
| **rho_R2G** | Spearman correlation between region activity and gene expression. |
| **importance_x_rho** | Combined metric: `importance_R2G × rho_R2G`. |
| **importance_x_abs_rho** | Combined metric with absolute correlation: `importance_R2G × |rho_R2G|`. |
| **TF** | Transcription factor predicted to regulate the gene via this region. |
| **is_extended** | Indicates whether the eRegulon includes extended regions (`True/False`). |
| **eRegulon_name** | Name of the TF-centered eRegulon module. |
| **Gene_signature_name** | Name of the **gene signature** subset linked to the TF. |
| **Region_signature_name** | Name of the **region signature** subset linked to the TF. |
| **importance_TF2G** | Importance of the **TF → Gene** regulation. |
| **regulation** | TF effect on gene: `1 = activating`, `-1 = repressing`. |
| **rho_TF2G** | Spearman correlation between TF activity/expression and target gene. |
| **triplet_rank** | Rank of the TF-Region-Gene triplet within the eRegulon (lower = higher priority). |

---

## How to Interpret a Row

| Region | Gene | TF | R2G importance | R2G rho | TF2G importance | Regulation | Triplet rank |
|--------|------|----|----------------|---------|----------------|-----------|--------------|
| chr2:71980454-71981615 | Rapgef4 | Alx1 | 0.067 | 0.111 | 3.21 | 1 (activating) | 10464 |

**Explanation:**  
- Region `chr2:71980454-71981615` is likely regulating **Rapgef4**.  
- **Alx1** is the TF controlling this regulation.  
- Both the region-to-gene and TF-to-gene links have moderate importance and correlation.  
- Triplet rank indicates relative strength within the eRegulon.

---

### Key Takeaways

1. Each row links **one region, one TF, and one target gene**.  
2. Metrics (`importance`, `rho`) help prioritize **strong regulatory relationships**.  
3. `is_extended` and signature columns provide context on **which regions/genes are considered core vs extended**.  
4. Triplets with **higher importance and lower rank** are likely the most biologically relevant.

---
## UMAP of SCENIC+ Regulon Activity by Cell Type

![UMAP of regulon activity by cell type](SnakemakeFast/workflow/figures/umap_cell_type.png)

*Figure:* UMAP projection of single cells based on **SCENIC+ regulon activity** rather than raw gene expression.  
Each point represents a single cell, and colors indicate cell types (e.g., AC, BC, Cones, MG, MGPC, Microglia, Rod).

---

### 🧠 What this plot shows
This is a **UMAP (Uniform Manifold Approximation and Projection)** built from the **regulon activity matrix** generated by **SCENIC+**.  
Instead of clustering cells based on expression of thousands of individual genes, this visualization captures **transcription factor (TF) activity patterns** inferred from gene regulatory networks.

In other words:
- Each cell’s position in this UMAP reflects **how similar its active TF regulons** are to those in other cells.  
- Cells with **similar regulatory programs** cluster together, forming distinct groups that correspond to biological cell types.


### 🔍 How to interpret
- The **tight clusters** indicate that each cell type has a **distinct regulon activity profile**, confirming that transcription factor programs strongly define cell identity.  
- **Proximity between clusters** reflects **regulatory similarity** — for example, Müller glia (MG) and Müller glia progenitor cells (MGPC) often lie near each other, consistent with their shared lineage and transcriptional programs.  
- The structure of the UMAP can reveal **transitions or relationships** between cell states — e.g., progenitor to mature cell trajectories.


#### UMAP of SCENIC+ eRegulon Activity: Top Features

This UMAP visualization represents the activity of the top regulatory features (eRegulons) across cells in our dataset, as inferred from SCENIC+.  

- Each point corresponds to a single cell.  
- Cells are colored based on the activity of selected top features: either the top 2 features per cell type from the `rss` ranking, or the top 10 most variable features if `rss` was not available.  
- Clustering patterns reflect similarities in regulatory activity rather than gene expression alone, highlighting potential cell type–specific regulatory programs.  
- This graph helps identify which eRegulons are most influential in driving cell type identity and can guide downstream analyses such as marker identification or regulatory network inference.


![UMAP of Top Features](SnakemakeFast/workflow/figures/umap_top_features.png)


#### SCENIC+ eRegulon Specificity Score (RSS) Plot

![SCENIC+ eRegulon Specificity Score (RSS) Plot](SnakemakeFast/workflow/RSS_plot.png)

This figure shows the **SCENIC+ eRegulon Specificity Score (RSS)** for each cell type.  
Each panel corresponds to one cell type or cluster (e.g., AC, BC, MG, MGPC, etc.), and each **dot** represents a **regulon** — a transcription factor (TF) and its predicted target genes.

### 🧠 What the axes mean
- **x-axis (eRegulon rank):** Regulons are ranked from most to least specific for that cell type.  
- **y-axis (eRegulon Specificity Score, RSS):** The higher the value, the more specifically that regulon is active in this cell type compared to others.  
  High RSS means the regulon is *enriched* and likely important for defining that cell type’s identity.

### 🔴 Red-labeled points
The red labels indicate the **top enriched regulons** — the TFs whose target genes are most specific and active in that cell type.  
Each label includes:
- **TF name** (e.g., *Emx1*, *Neurod2*, *Otx2*)  
- Whether it’s **direct** (+/+) or **indirect** (+/–, –/–), based on motif evidence  
- The **number of target genes** (e.g., 76g = 76 genes)

Example:  
`Emx1_direct_+/+_(76g)` → *Emx1* regulon, direct motif evidence, 76 target genes.

### 📊 How to interpret
The **higher** a regulon appears, the more **specific** and **enriched** it is in that cell type.  
These top TFs are likely the *key transcriptional regulators* defining each cell population.


## References

- [Pycistopic: Human Cerebellum Notebook](https://pycistopic.readthedocs.io/en/latest/notebooks/human_cerebellum.html#Getting-pseudobulk-profiles-from-cell-annotations)
- [SCENIC+: Official Documentation](https://scenicplus.readthedocs.io/en/latest/index.html)



