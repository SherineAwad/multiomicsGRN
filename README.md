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
![Dotplot](figures/clustered_mNeurog2_Dotplot.png)

### Marker Genes UMAP

<img src="figures/umap_clustered_mNeurog2_Abca8a.png?v=3" alt="Abca8a" width="33%"><img src="figures/umap_clustered_mNeurog2_Ccr2.png?v=3" alt="Ccr2" width="33%"><img src="figures/umap_clustered_mNeurog2_Hes1.png?v=3" alt="Hes1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Notch1.png?v=3" alt="Notch1" width="33%"><img src="figures/umap_clustered_mNeurog2_Rlbp1.png?v=3" alt="Rlbp1" width="33%"><img src="figures/umap_clustered_mNeurog2_Acta2.png?v=3" alt="Acta2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Chat.png?v=3" alt="Chat" width="33%"><img src="figures/umap_clustered_mNeurog2_Hes5.png?v=3" alt="Hes5" width="33%"><img src="figures/umap_clustered_mNeurog2_Nrl.png?v=3" alt="Nrl" width="33%">
<img src="figures/umap_clustered_mNeurog2_Rpe65.png?v=3" alt="Rpe65" width="33%"><img src="figures/umap_clustered_mNeurog2_Apoe.png?v=3" alt="Apoe" width="33%">
<img src="figures/umap_clustered_mNeurog2_Insm1.png?v=3" alt="Insm1" width="33%"><img src="figures/umap_clustered_mNeurog2_Olig2.png?v=3" alt="Olig2" width="33%"><img src="figures/umap_clustered_mNeurog2_Sebox.png?v=3" alt="Sebox" width="33%">
<img src="figures/umap_clustered_mNeurog2_Aqp4.png?v=3" alt="Aqp4" width="33%"><img src="figures/umap_clustered_mNeurog2_Csf1r.png?v=3" alt="Csf1r" width="33%"><img src="figures/umap_clustered_mNeurog2_Isl1.png?v=3" alt="Isl1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Otx2.png?v=3" alt="Otx2" width="33%"><img src="figures/umap_clustered_mNeurog2_Slc17a7.png?v=3" alt="Slc17a7" width="33%"><img src="figures/umap_clustered_mNeurog2_Arr3.png?v=3" alt="Arr3" width="33%">
<img src="figures/umap_clustered_mNeurog2_Elavl3.png?v=3" alt="Elavl3" width="33%"><img src="figures/umap_clustered_mNeurog2_Kcnj8.png?v=3" alt="Kcnj8" width="33%"><img src="figures/umap_clustered_mNeurog2_Pax2.png?v=3" alt="Pax2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Slc1a3.png?v=3" alt="Slc1a3" width="33%"><img src="figures/umap_clustered_mNeurog2_Ascl1.png?v=3" alt="Ascl1" width="33%"><img src="figures/umap_clustered_mNeurog2_Elavl4.png?v=3" alt="Elavl4" width="33%">
<img src="figures/umap_clustered_mNeurog2_Lhx1.png?v=3" alt="Lhx1" width="33%"><img src="figures/umap_clustered_mNeurog2_Pax6.png?v=3" alt="Pax6" width="33%"><img src="figures/umap_clustered_mNeurog2_Slc6a9.png?v=3" alt="Slc6a9" width="33%">
<img src="figures/umap_clustered_mNeurog2_Atoh7.png?v=3" alt="Atoh7" width="33%"><img src="figures/umap_clustered_mNeurog2_Emx1.png?v=3" alt="Emx1" width="33%"><img src="figures/umap_clustered_mNeurog2_Lhx2.png?v=3" alt="Lhx2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Pou4f2.png?v=3" alt="Pou4f2" width="33%"><img src="figures/umap_clustered_mNeurog2_Sox11.png?v=3" alt="Sox11" width="33%"><img src="figures/umap_clustered_mNeurog2_Bsn.png?v=3" alt="Bsn" width="33%">
<img src="figures/umap_clustered_mNeurog2_Foxn4.png?v=3" alt="Foxn4" width="33%"><img src="figures/umap_clustered_mNeurog2_Lhx4.png?v=3" alt="Lhx4" width="33%"><img src="figures/umap_clustered_mNeurog2_Prdm1.png?v=3" alt="Prdm1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Sox9.png?v=3" alt="Sox9" width="33%"><img src="figures/umap_clustered_mNeurog2_Cabp5.png?v=3" alt="Cabp5" width="33%"><img src="figures/umap_clustered_mNeurog2_Gad1.png?v=3" alt="Gad1" width="33%">
<img src="figures/umap_clustered_mNeurog2_Malat1.png?v=3" alt="Malat1" width="33%"><img src="figures/umap_clustered_mNeurog2_Prdx6.png?v=3" alt="Prdx6" width="33%"><img src="figures/umap_clustered_mNeurog2_Tfap2a.png?v=3" alt="Tfap2a" width="33%">
<img src="figures/umap_clustered_mNeurog2_Calb1.png?v=3" alt="Calb1" width="33%"><img src="figures/umap_clustered_mNeurog2_Gfap.png?v=3" alt="Gfap" width="33%"><img src="figures/umap_clustered_mNeurog2_mt-Atp6.png?v=3" alt="mt-Atp6" width="33%">
<img src="figures/umap_clustered_mNeurog2_Rbfox3.png?v=3" alt="Rbfox3" width="33%"><img src="figures/umap_clustered_mNeurog2_Tie1.png?v=3" alt="Tie1" width="33%"><img src="figures/umap_clustered_mNeurog2_Calb2.png?v=3" alt="Calb2" width="33%">
<img src="figures/umap_clustered_mNeurog2_Glul.png?v=3" alt="Glul" width="33%"><img src="figures/umap_clustered_mNeurog2_Neurog2.png?v=3" alt="Neurog2" width="33%"><img src="figures/umap_clustered_mNeurog2_Rho.png?v=3" alt="Rho" width="33%">
<img src="figures/umap_clustered_mNeurog2_Vim.png?v=3" alt="Vim" width="33%">

### Annotations 

![Annotations ON](figures/annotated_clustered_mNeurog2_annotationsON.png)

![Annotations](figures/annotated_clustered_mNeurog2_annotations.png)


# Part 2: Analyzing scATAC Data  

![ATAC UMAP](ATAC_samples.png)

### Clustering 

![ATAC CLUSTERS](ATAC_clusters.png)


### ATAC Markers UMAP 

<img src="ATAC_FeaturePlot_Abca8a.png" width="33%"><img src="ATAC_FeaturePlot_Acta2.png" width="33%"><img src="ATAC_FeaturePlot_Apoe.png" width="33%">
<img src="ATAC_FeaturePlot_Aqp4.png" width="33%"><img src="ATAC_FeaturePlot_Arr3.png" width="33%"><img src="ATAC_FeaturePlot_Ascl1.png" width="33%">
<img src="ATAC_FeaturePlot_Atoh7.png" width="33%"><img src="ATAC_FeaturePlot_Bhlhe23.png" width="33%"><img src="ATAC_FeaturePlot_Bsn.png" width="33%">
<img src="ATAC_FeaturePlot_Cabp5.png" width="33%"><img src="ATAC_FeaturePlot_Calb1.png" width="33%"><img src="ATAC_FeaturePlot_Calb2.png" width="33%">
<img src="ATAC_FeaturePlot_Cbln4.png" width="33%"><img src="ATAC_FeaturePlot_Ccr2.png" width="33%"><img src="ATAC_FeaturePlot_Chat.png" width="33%">
<img src="ATAC_FeaturePlot_Crx.png" width="33%"><img src="ATAC_FeaturePlot_Csf1r.png" width="33%"><img src="ATAC_FeaturePlot_Cx3cr1.png" width="33%">
<img src="ATAC_FeaturePlot_Ebf3.png" width="33%"><img src="ATAC_FeaturePlot_Elavl3.png" width="33%"><img src="ATAC_FeaturePlot_Elavl4.png" width="33%">
<img src="ATAC_FeaturePlot_Emx1.png" width="33%"><img src="ATAC_FeaturePlot_Fgf15.png" width="33%"><img src="ATAC_FeaturePlot_Foxn4.png" width="33%">
<img src="ATAC_FeaturePlot_Gad1.png" width="33%"><img src="ATAC_FeaturePlot_Gad2.png" width="33%"><img src="ATAC_FeaturePlot_Gfap.png" width="33%">
<img src="ATAC_FeaturePlot_Gli1.png" width="33%"><img src="ATAC_FeaturePlot_Glul.png" width="33%"><img src="ATAC_FeaturePlot_Gnat2.png" width="33%">
<img src="ATAC_FeaturePlot_Grm6.png" width="33%"><img src="ATAC_FeaturePlot_Guca1b.png" width="33%"><img src="ATAC_FeaturePlot_Hes1.png" width="33%">
<img src="ATAC_FeaturePlot_Hes5.png" width="33%"><img src="ATAC_FeaturePlot_Igf2.png" width="33%"><img src="ATAC_FeaturePlot_Insm1.png" width="33%">
<img src="ATAC_FeaturePlot_Isl1.png" width="33%"><img src="ATAC_FeaturePlot_Isl2.png" width="33%"><img src="ATAC_FeaturePlot_Kcnj8.png" width="33%">
<img src="ATAC_FeaturePlot_Lhx1.png" width="33%"><img src="ATAC_FeaturePlot_Lhx2.png" width="33%"><img src="ATAC_FeaturePlot_Lhx4.png" width="33%">
<img src="ATAC_FeaturePlot_mt-Atp6.png" width="33%"><img src="ATAC_FeaturePlot_Nefl.png" width="33%"><img src="ATAC_FeaturePlot_Nefm.png" width="33%">
<img src="ATAC_FeaturePlot_Neurog2.png" width="33%"><img src="ATAC_FeaturePlot_Notch1.png" width="33%"><img src="ATAC_FeaturePlot_Nr2e3.png" width="33%">
<img src="ATAC_FeaturePlot_Olig2.png" width="33%"><img src="ATAC_FeaturePlot_Onecut1.png" width="33%"><img src="ATAC_FeaturePlot_Onecut2.png" width="33%">
<img src="ATAC_FeaturePlot_Opn1mw.png" width="33%"><img src="ATAC_FeaturePlot_Opn1sw.png" width="33%"><img src="ATAC_FeaturePlot_Otx2.png" width="33%">
<img src="ATAC_FeaturePlot_Pax2.png" width="33%"><img src="ATAC_FeaturePlot_Pax6.png" width="33%"><img src="ATAC_FeaturePlot_Pcp4.png" width="33%">
<img src="ATAC_FeaturePlot_Pdgfra.png" width="33%"><img src="ATAC_FeaturePlot_Pecam1.png" width="33%"><img src="ATAC_FeaturePlot_Pou4f1.png" width="33%">
<img src="ATAC_FeaturePlot_Pou4f2.png" width="33%"><img src="ATAC_FeaturePlot_Pou4f3.png" width="33%"><img src="ATAC_FeaturePlot_Prdm1.png" width="33%">
<img src="ATAC_FeaturePlot_Prdx6.png" width="33%"><img src="ATAC_FeaturePlot_Prkca.png" width="33%"><img src="ATAC_FeaturePlot_Prox1.png" width="33%">
<img src="ATAC_FeaturePlot_Ptprc.png" width="33%"><img src="ATAC_FeaturePlot_Rbfox3.png" width="33%"><img src="ATAC_FeaturePlot_Rbpms.png" width="33%">
<img src="ATAC_FeaturePlot_Rho.png" width="33%"><img src="ATAC_FeaturePlot_Rlbp1.png" width="33%"><img src="ATAC_FeaturePlot_Rom1.png" width="33%">
<img src="ATAC_FeaturePlot_Sall1.png" width="33%"><img src="ATAC_FeaturePlot_Sebox.png" width="33%"><img src="ATAC_FeaturePlot_Sfrp2.png" width="33%">
<img src="ATAC_FeaturePlot_Slc17a7.png" width="33%"><img src="ATAC_FeaturePlot_Slc18a3.png" width="33%"><img src="ATAC_FeaturePlot_Slc1a3.png" width="33%">
<img src="ATAC_FeaturePlot_Slc6a9.png" width="33%"><img src="ATAC_FeaturePlot_Sncg.png" width="33%"><img src="ATAC_FeaturePlot_Sox11.png" width="33%">
<img src="ATAC_FeaturePlot_Sox2.png" width="33%"><img src="ATAC_FeaturePlot_Sox9.png" width="33%"><img src="ATAC_FeaturePlot_Tfap2a.png" width="33%">
<img src="ATAC_FeaturePlot_Tfap2b.png" width="33%"><img src="ATAC_FeaturePlot_Thrb.png" width="33%"><img src="ATAC_FeaturePlot_Thy1.png" width="33%">
<img src="ATAC_FeaturePlot_Tie1.png" width="33%"><img src="ATAC_FeaturePlot_Trpm1.png" width="33%"><img src="ATAC_FeaturePlot_Vim.png" width="33%">
<img src="ATAC_FeaturePlot_Vsx1.png" width="33%"><img src="ATAC_FeaturePlot_Vsx2.png" width="33%">



#### Barcode-to-Cluster Mapping

The following file format (`CSV`) shows how barcodes are mapped to clusters:
After annotation, we can replace the id with the cluster 


#### Header of the barcode to cluster mapping 
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

