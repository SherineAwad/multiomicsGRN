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

<img src="figures/umap_clustered_mNeurog2_Abca8a.png"> <img src="figures/umap_clustered_mNeurog2_Ccr2.png"> <img src="figures/umap_clustered_mNeurog2_Hes1.png">  
<img src="figures/umap_clustered_mNeurog2_Notch1.png"> <img src="figures/umap_clustered_mNeurog2_Rlbp1.png"> <img src="figures/umap_clustered_mNeurog2_Acta2.png">  
<img src="figures/umap_clustered_mNeurog2_Chat.png"> <img src="figures/umap_clustered_mNeurog2_Hes5.png"> <img src="figures/umap_clustered_mNeurog2_Nrl.png">  
<img src="figures/umap_clustered_mNeurog2_Rpe65.png"> <img src="figures/umap_clustered_mNeurog2_Apoe.png"> <img src="figures/umap_clustered_mNeurog2_Clusters.png">  
<img src="figures/umap_clustered_mNeurog2_Insm1.png"> <img src="figures/umap_clustered_mNeurog2_Olig2.png"> <img src="figures/umap_clustered_mNeurog2_Sebox.png">  
<img src="figures/umap_clustered_mNeurog2_Aqp4.png"> <img src="figures/umap_clustered_mNeurog2_Csf1r.png"> <img src="figures/umap_clustered_mNeurog2_Isl1.png">  
<img src="figures/umap_clustered_mNeurog2_Otx2.png"> <img src="figures/umap_clustered_mNeurog2_Slc17a7.png"> <img src="figures/umap_clustered_mNeurog2_Arr3.png">  
<img src="figures/umap_clustered_mNeurog2_Elavl3.png"> <img src="figures/umap_clustered_mNeurog2_Kcnj8.png"> <img src="figures/umap_clustered_mNeurog2_Pax2.png">  
<img src="figures/umap_clustered_mNeurog2_Slc1a3.png"> <img src="figures/umap_clustered_mNeurog2_Ascl1.png"> <img src="figures/umap_clustered_mNeurog2_Elavl4.png">  
<img src="figures/umap_clustered_mNeurog2_Lhx1.png"> <img src="figures/umap_clustered_mNeurog2_Pax6.png"> <img src="figures/umap_clustered_mNeurog2_Slc6a9.png">  
<img src="figures/umap_clustered_mNeurog2_Atoh7.png"> <img src="figures/umap_clustered_mNeurog2_Emx1.png"> <img src="figures/umap_clustered_mNeurog2_Lhx2.png">  
<img src="figures/umap_clustered_mNeurog2_Pou4f2.png"> <img src="figures/umap_clustered_mNeurog2_Sox11.png"> <img src="figures/umap_clustered_mNeurog2_Bsn.png">  
<img src="figures/umap_clustered_mNeurog2_Foxn4.png"> <img src="figures/umap_clustered_mNeurog2_Lhx4.png"> <img src="figures/umap_clustered_mNeurog2_Prdm1.png">  
<img src="figures/umap_clustered_mNeurog2_Sox9.png"> <img src="figures/umap_clustered_mNeurog2_Cabp5.png"> <img src="figures/umap_clustered_mNeurog2_Gad1.png">  
<img src="figures/umap_clustered_mNeurog2_Malat1.png"> <img src="figures/umap_clustered_mNeurog2_Prdx6.png"> <img src="figures/umap_clustered_mNeurog2_Tfap2a.png">  
<img src="figures/umap_clustered_mNeurog2_Calb1.png"> <img src="figures/umap_clustered_mNeurog2_Gfap.png"> <img src="figures/umap_clustered_mNeurog2_mt-Atp6.png">  
<img src="figures/umap_clustered_mNeurog2_Rbfox3.png"> <img src="figures/umap_clustered_mNeurog2_Tie1.png"> <img src="figures/umap_clustered_mNeurog2_Calb2.png">  
<img src="figures/umap_clustered_mNeurog2_Glul.png"> <img src="figures/umap_clustered_mNeurog2_Neurog2.png"> <img src="figures/umap_clustered_mNeurog2_Rho.png">  
<img src="figures/umap_clustered_mNeurog2_Vim.png">




