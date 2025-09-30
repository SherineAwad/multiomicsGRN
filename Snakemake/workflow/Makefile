
chromsizes.tsv:
	samtools faidx /nfs/turbo/umms-thahoang/sherine/REFERENCES/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa
	cut -f1,2 /nfs/turbo/umms-thahoang/sherine/REFERENCES/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa.fai > chromsizes.tsv



GENOME_GTF = /nfs/turbo/umms-thahoang/sherine/REFERENCES/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf
CHROMSIZES = /nfs/turbo/umms-thahoang/sherine/scenicplus_neurog2/chromsizes.tsv
GENOME_ANNOTATION = genome_annotation.tsv

GENOME_ANNOTATION = genome_annotation.tsv
GTF = /nfs/turbo/umms-thahoang/sherine/REFERENCES/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf

genome_annotation.tsv:
	bash make_annot_file.sh


