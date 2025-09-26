

mm10-blacklist.v2.bed:
	wget https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz
	gunzip ENCFF547MET.bed.gz
	mv ENCFF547MET.bed mm10-blacklist.v2.bed


## Set up snakemake scenicplus init_snakemake --out_dir scplus_pipeline
motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl:
	wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl

##https://resources.aertslab.org/cistarget/
##https://resources.aertslab.org/cistarget/databases/


mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather:
	wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather


mm10_screen_v10_clust.regions_vs_motifs.rankings.feather:
	wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather
