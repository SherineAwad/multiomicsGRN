# Multi-omics GRN Construction with SCENIC+  

We are using **multi-omics data for Neurog2** to construct a gene regulatory network (GRN) with **SCENIC+**.  

# Part 1: Analyzing scRNA-seq Data  

![Before Filtering QC metrics](figures/violin_QC.png)
Violin plots displaying quality control metrics such as number of genes detected per cell, total counts, and percentage of mitochondrial gene expression.

### Filtering Criteria

Quality filtering was applied to remove low-quality cells and potential doublets. Cells were retained only if they met all the following criteria:

- Number of genes detected per cell between **800 and 8000**
- Total counts per cell between **1200 and 30000**
- Percentage of mitochondrial gene counts less than **25%**

This filtering step ensures removal of dead or dying cells and technical artifacts to improve downstream analysis quality.


### Additional Analysis Figure
![After Filtering QC metrics](figures/violin_AfterQC.png)

### UMAP and Clusters UMAP


![Samples UMAP](figures/umap_clustered_mNeurog2_Samples.png)


![Clusters UMAP](figures/umap_clustered_mNeurog2_Clusters.png)

### Marker Genes UMAP 

### Marker Genes UMAP

<img src="figures/umap_clustered_mNeurog2_Abca8a.png" alt="Abca8a" width="33%"><img src="figures/umap_clustered_mNeurog2_Ccr2.png" alt="Ccr2" width="33%"><img src="figures/umap_clustered_mNeurog2_Hes1.png" alt="Hes1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Notch1.png" alt="Notch1" width="33%"><img src="figures/umap_clustered_mNeurog2_Rlbp1.png" alt="Rlbp1" width="33%"><img src="figures/umap_clustered_mNeurog2_Acta2.png" alt="Acta2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Chat.png" alt="Chat" width="33%"><img src="figures/umap_clustered_mNeurog2_Hes5.png" alt="Hes5" width="33%"><img src="figures/umap_clustered_mNeurog2_Nrl.png" alt="Nrl" width="33%">
<img src="figures/umap_clustered_mNeurog2_Rpe65.png" alt="Rpe65" width="33%"><img src="figures/umap_clustered_mNeurog2_Apoe.png" alt="Apoe" width="33%"><img src="figures/umap_clustered_mNeurog2_Clusters.png" alt="Clusters" width="33%">
<img src="figures/umap_clustered_mNeurog2_Insm1.png" alt="Insm1" width="33%"><img src="figures/umap_clustered_mNeurog2_Olig2.png" alt="Olig2" width="33%"><img src="figures/umap_clustered_mNeurog2_Sebox.png" alt="Sebox" width="33%">
<img src="figures/umap_clustered_mNeurog2_Aqp4.png" alt="Aqp4" width="33%"><img src="figures/umap_clustered_mNeurog2_Csf1r.png" alt="Csf1r" width="33%"><img src="figures/umap_clustered_mNeurog2_Isl1.png" alt="Isl1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Otx2.png" alt="Otx2" width="33%"><img src="figures/umap_clustered_mNeurog2_Slc17a7.png" alt="Slc17a7" width="33%"><img src="figures/umap_clustered_mNeurog2_Arr3.png" alt="Arr3" width="33%">
<img src="figures/umap_clustered_mNeurog2_Elavl3.png" alt="Elavl3" width="33%"><img src="figures/umap_clustered_mNeurog2_Kcnj8.png" alt="Kcnj8" width="33%"><img src="figures/umap_clustered_mNeurog2_Pax2.png" alt="Pax2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Slc1a3.png" alt="Slc1a3" width="33%"><img src="figures/umap_clustered_mNeurog2_Ascl1.png" alt="Ascl1" width="33%"><img src="figures/umap_clustered_mNeurog2_Elavl4.png" alt="Elavl4" width="33%">
<img src="figures/umap_clustered_mNeurog2_Lhx1.png" alt="Lhx1" width="33%"><img src="figures/umap_clustered_mNeurog2_Pax6.png" alt="Pax6" width="33%"><img src="figures/umap_clustered_mNeurog2_Slc6a9.png" alt="Slc6a9" width="33%">
<img src="figures/umap_clustered_mNeurog2_Atoh7.png" alt="Atoh7" width="33%"><img src="figures/umap_clustered_mNeurog2_Emx1.png" alt="Emx1" width="33%"><img src="figures/umap_clustered_mNeurog2_Lhx2.png" alt="Lhx2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Pou4f2.png" alt="Pou4f2" width="33%"><img src="figures/umap_clustered_mNeurog2_Sox11.png" alt="Sox11" width="33%"><img src="figures/umap_clustered_mNeurog2_Bsn.png" alt="Bsn" width="33%">
<img src="figures/umap_clustered_mNeurog2_Foxn4.png" alt="Foxn4" width="33%"><img src="figures/umap_clustered_mNeurog2_Lhx4.png" alt="Lhx4" width="33%"><img src="figures/umap_clustered_mNeurog2_Prdm1.png" alt="Prdm1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Sox9.png" alt="Sox9" width="33%"><img src="figures/umap_clustered_mNeurog2_Cabp5.png" alt="Cabp5" width="33%"><img src="figures/umap_clustered_mNeurog2_Gad1.png" alt="Gad1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Malat1.png" alt="Malat1" width="33%"><img src="figures/umap_clustered_mNeurog2_Prdx6.png" alt="Prdx6" width="33%"><img src="figures/umap_clustered_mNeurog2_Tfap2a.png" alt="Tfap2a" width="33%">
<img src="figures/umap_clustered_mNeurog2_Calb1.png" alt="Calb1" width="33%"><img src="figures/umap_clustered_mNeurog2_Gfap.png" alt="Gfap" width="33%"><img src="figures/umap_clustered_mNeurog2_mt-Atp6.png" alt="mt-Atp6" width="33%">
<img src="figures/umap_clustered_mNeurog2_Rbfox3.png" alt="Rbfox3" width="33%"><img src="figures/umap_clustered_mNeurog2_Tie1.png" alt="Tie1" width="33%"><img src="figures/umap_clustered_mNeurog2_Calb2.png" alt="Calb2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Glul.png" alt="Glul" width="33%"><img src="figures/umap_clustered_mNeurog2_Neurog2.png" alt="Neurog2" width="33%"><img src="figures/umap_clustered_mNeurog2_Rho.png" alt="Rho" width="33%">
<img src="figures/umap_clustered_mNeurog2_Vim.png" alt="Vim" width="33%">





# Part 2: Analyzing scATAC Data  

![ATAC UMAP](atac_umap.png)

###GOING ON NEED MORE PLOTS TO ANNOTATE

## Barcode-to-Cluster Mapping

The following file format (`CSV`) shows how barcodes are mapped to clusters:
After annotation, we can replace the id with the cluster 


## Header of the barcode to cluster mapping 
```csv
"barcode","cluster"
"Control_CGCCTCATCCTAAGGT-1","4"
"Control_CCAGCCTGTTTGGGCG-1","0"
"Control_TTGACTAAGTCATTTC-1","2"
"Control_GACTATTCAAGGTCCT-1","0"
"Control_GCTAGCTCACAGGAAT-1","10"
"Control_AGCATCCCACCATATG-1","0"
"Control_AGGTTAGAGCGATACT-1","11"
"Control_TCTCACCAGCTTCTCA-1","2"
"Control_GGCTCAATCCTAGTTT-1","0"
"Control_TTGCATTTCTCTAGCC-1","10"
"Control_CACCTCAGTGTTTGAG-1","8"
"Control_AGTAGGATCACTCAAA-1","0"
"Control_CCTAAAGGTAAGGTTT-1","4"
"Control_GTGCTTACAGCCAGAA-1","0"
"Control_CTCTTGATCGAGGTGG-1","10"
"Control_CCTACTTCACATTGCA-1","0"
"Control_GCCTACTTCTGTAATG-1","4"
"Control_AATCCGTAGTATCGCG-1","0"
"Control_GTGGATGCAAGACTCC-1","10"

