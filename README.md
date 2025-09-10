# Multi-omics GRN Construction with SCENIC+  

We are using **multi-omics data for Neurog2** to construct a gene regulatory network (GRN) with **SCENIC+**.  

# Part 1: Analyzing scRNA-seq Data  

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
![Dotplot](figures/Dotplot.png) 

### Marker Genes UMAP

<img src="figures/umap_clustered_mNeurog2_Abca8a.png?v=2" alt="Abca8a" width="33%"><img src="figures/umap_clustered_mNeurog2_Ccr2.png?v=2" alt="Ccr2" width="33%"><img src="figures/umap_clustered_mNeurog2_Hes1.png?v=2" alt="Hes1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Notch1.png?v=2" alt="Notch1" width="33%"><img src="figures/umap_clustered_mNeurog2_Rlbp1.png?v=2" alt="Rlbp1" width="33%"><img src="figures/umap_clustered_mNeurog2_Acta2.png?v=2" alt="Acta2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Chat.png?v=2" alt="Chat" width="33%"><img src="figures/umap_clustered_mNeurog2_Hes5.png?v=2" alt="Hes5" width="33%"><img src="figures/umap_clustered_mNeurog2_Nrl.png?v=2" alt="Nrl" width="33%">
<img src="figures/umap_clustered_mNeurog2_Rpe65.png?v=2" alt="Rpe65" width="33%"><img src="figures/umap_clustered_mNeurog2_Apoe.png?v=2" alt="Apoe" width="33%">
<img src="figures/umap_clustered_mNeurog2_Insm1.png?v=2" alt="Insm1" width="33%"><img src="figures/umap_clustered_mNeurog2_Olig2.png?v=2" alt="Olig2" width="33%"><img src="figures/umap_clustered_mNeurog2_Sebox.png?v=2" alt="Sebox" width="33%">
<img src="figures/umap_clustered_mNeurog2_Aqp4.png?v=2" alt="Aqp4" width="33%"><img src="figures/umap_clustered_mNeurog2_Csf1r.png?v=2" alt="Csf1r" width="33%"><img src="figures/umap_clustered_mNeurog2_Isl1.png?v=2" alt="Isl1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Otx2.png?v=2" alt="Otx2" width="33%"><img src="figures/umap_clustered_mNeurog2_Slc17a7.png?v=2" alt="Slc17a7" width="33%"><img src="figures/umap_clustered_mNeurog2_Arr3.png?v=2" alt="Arr3" width="33%">
<img src="figures/umap_clustered_mNeurog2_Elavl3.png?v=2" alt="Elavl3" width="33%"><img src="figures/umap_clustered_mNeurog2_Kcnj8.png?v=2" alt="Kcnj8" width="33%"><img src="figures/umap_clustered_mNeurog2_Pax2.png?v=2" alt="Pax2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Slc1a3.png?v=2" alt="Slc1a3" width="33%"><img src="figures/umap_clustered_mNeurog2_Ascl1.png?v=2" alt="Ascl1" width="33%"><img src="figures/umap_clustered_mNeurog2_Elavl4.png?v=2" alt="Elavl4" width="33%">
<img src="figures/umap_clustered_mNeurog2_Lhx1.png?v=2" alt="Lhx1" width="33%"><img src="figures/umap_clustered_mNeurog2_Pax6.png?v=2" alt="Pax6" width="33%"><img src="figures/umap_clustered_mNeurog2_Slc6a9.png?v=2" alt="Slc6a9" width="33%">
<img src="figures/umap_clustered_mNeurog2_Atoh7.png?v=2" alt="Atoh7" width="33%"><img src="figures/umap_clustered_mNeurog2_Emx1.png?v=2" alt="Emx1" width="33%"><img src="figures/umap_clustered_mNeurog2_Lhx2.png?v=2" alt="Lhx2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Pou4f2.png?v=2" alt="Pou4f2" width="33%"><img src="figures/umap_clustered_mNeurog2_Sox11.png?v=2" alt="Sox11" width="33%"><img src="figures/umap_clustered_mNeurog2_Bsn.png?v=2" alt="Bsn" width="33%">
<img src="figures/umap_clustered_mNeurog2_Foxn4.png?v=2" alt="Foxn4" width="33%"><img src="figures/umap_clustered_mNeurog2_Lhx4.png?v=2" alt="Lhx4" width="33%"><img src="figures/umap_clustered_mNeurog2_Prdm1.png?v=2" alt="Prdm1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Sox9.png?v=2" alt="Sox9" width="33%"><img src="figures/umap_clustered_mNeurog2_Cabp5.png?v=2" alt="Cabp5" width="33%"><img src="figures/umap_clustered_mNeurog2_Gad1.png?v=2" alt="Gad1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Malat1.png?v=2" alt="Malat1" width="33%"><img src="figures/umap_clustered_mNeurog2_Prdx6.png?v=2" alt="Prdx6" width="33%"><img src="figures/umap_clustered_mNeurog2_Tfap2a.png?v=2" alt="Tfap2a" width="33%">
<img src="figures/umap_clustered_mNeurog2_Calb1.png?v=2" alt="Calb1" width="33%"><img src="figures/umap_clustered_mNeurog2_Gfap.png?v=2" alt="Gfap" width="33%"><img src="figures/umap_clustered_mNeurog2_mt-Atp6.png?v=2" alt="mt-Atp6" width="33%">
<img src="figures/umap_clustered_mNeurog2_Rbfox3.png?v=2" alt="Rbfox3" width="33%"><img src="figures/umap_clustered_mNeurog2_Tie1.png?v=2" alt="Tie1" width="33%"><img src="figures/umap_clustered_mNeurog2_Calb2.png?v=2" alt="Calb2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Glul.png?v=2" alt="Glul" width="33%"><img src="figures/umap_clustered_mNeurog2_Neurog2.png?v=2" alt="Neurog2" width="33%"><img src="figures/umap_clustered_mNeurog2_Rho.png?v=2" alt="Rho" width="33%">
<img src="figures/umap_clustered_mNeurog2_Vim.png?v=2" alt="Vim" width="33%">



# Part 2: Analyzing scATAC Data  

![ATAC UMAP](ATAC_samples.png)

###Clustering 

![ATAC CLUSTERS](ATAC_clusters.png)


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

