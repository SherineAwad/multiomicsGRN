import os
configfile: "config.yaml"

RESULTS = config["paths"]["results_dir"]

rule all:
    input:
        f"{RESULTS}/cistopic_objects_mm10.pkl",
        f"{RESULTS}/merged_cistopic.pkl",
        f"{RESULTS}/merged_with_meta.pkl",
        f"{RESULTS}/MALLET/merged_cistopic_with_models.pkl",
        f"{RESULTS}/umap_clusters/cistopic_obj_clustered.pkl",
        f"{RESULTS}/topics/cistopic_obj_binarized.pkl",
        f"{RESULTS}/DAR_results/cistopic_obj_with_DARs.pkl"

rule pycistopic_pseudobulk:
    shell:
        "python pycistopic_pseudobulk.py"

rule run_macs2:
    shell:
        "python run_macs2.py"

rule consensus_peaks:
    input:
        "all_peaks_combined.bed"
    output:
        f"{RESULTS}/consensus_peak_calling/consensus_peaks.bed"
    shell:
        "python consensus_peaks.py {input} {output}"

rule get_tss:
    output:
        f"{RESULTS}/QC/tss_mm10.bed"
    params:
        name="mmusculus_gene_ensembl",
        ucsc="mm10"
    shell:
        """
        mkdir -p {RESULTS}/QC
        pycistopic tss get_tss --output {output} --name {params.name} --ucsc {params.ucsc}
        """

rule run_qc:
    input:
        th1=config["fragments"]["th1"],
        th2=config["fragments"]["th2"],
        tss_bed=f"{RESULTS}/QC/tss_mm10.bed"
    output:
        directory=f"{RESULTS}/QC"
    shell:
        "python run_qc.py --out_dir {RESULTS} "
        "--consensus_dir {RESULTS}/consensus_peak_calling "
        "--tss_bed {input.tss_bed} "
        "--th1_fragments {input.th1} "
        "--th2_fragments {input.th2}"

rule pycistopic_qc_commands:
    input:
        "pycistopic_qc_commands.txt"
    shell:
        "bash {input}"

rule collect_qc_barcodes:
    input:
        fragments_dict=f"{RESULTS}/fragments_dict.pkl"
    output:
        pickle=f"{RESULTS}/QC/qc_barcodes_thresholds.pkl"
    params:
        unique_fragments_threshold=config["qc"]["unique_fragments_threshold"],
        tss_enrichment_threshold=config["qc"]["tss_enrichment_threshold"],
        frip_threshold=config["qc"]["frip_threshold"]
    shell:
        "python collect_qc_barcodes.py "
        "--fragments_dict {input.fragments_dict} "
        "--qc_output_dir {RESULTS}/QC "
        "--output_pickle {output.pickle} "
        "--unique_fragments_threshold {params.unique_fragments_threshold} "
        "--tss_enrichment_threshold {params.tss_enrichment_threshold} "
        "--frip_threshold {params.frip_threshold}"

rule plot_pycistopicQC:
    input:
        fragments_dict=f"{RESULTS}/fragments_dict.pkl",
        qc_results_pickle=f"{RESULTS}/QC/qc_barcodes_thresholds.pkl"
    output:
        plots_dir=f"{RESULTS}/QC"
    shell:
        "python plot_pycistopicQC.py "
        "--fragments_dict {input.fragments_dict} "
        "--qc_output_dir {output.plots_dir} "
        "--plots_output_dir {output.plots_dir} "
        "--qc_results_pickle {input.qc_results_pickle} "
        "--barcode_plots_output_dir {output.plots_dir}"

rule create_cistopic_objects:
    input:
        fragments_dict=f"{RESULTS}/fragments_dict.pkl",
        qc_results_pickle=f"{RESULTS}/QC/qc_barcodes_thresholds.pkl",
        regions_bed=f"{RESULTS}/consensus_peak_calling/consensus/consensus_peaks.bed",
        blacklist=config["paths"]["blacklist_bed"]
    output:
        pickle=f"{RESULTS}/cistopic_objects_mm10.pkl"
    params:
        n_cpu=config["resources"]["pycistopic_n_cpu"]
    shell:
        "python create_cistopic_objects.py "
        "--fragments_dict {input.fragments_dict} "
        "--qc_results_pickle {input.qc_results_pickle} "
        "--regions_bed {input.regions_bed} "
        "--blacklist_bed {input.blacklist} "
        "--qc_output_dir {RESULTS}/QC "
        "--output_pickle {output.pickle} "
        "--n_cpu {params.n_cpu}"

rule merge_cistopic:
    input:
        f"{RESULTS}/cistopic_objects_mm10.pkl"
    output:
        f"{RESULTS}/merged_cistopic.pkl"
    shell:
        "python merge_cistopic.py {input} {output}"

rule add_scrna_metadata:
    input:
        cistopic_pickle=f"{RESULTS}/merged_cistopic.pkl",
        scrna_csv=config["paths"]["scrna_csv"]
    output:
        f"{RESULTS}/merged_with_meta.pkl"
    shell:
        "python add_scrna_metadata_fixed.py "
        "--cistopic_pickle {input.cistopic_pickle} "
        "--scrna_csv {input.scrna_csv} "
        "--output_pickle {output}"

rule run_mallet:
    input:
        cistopic_pickle=f"{RESULTS}/merged_with_meta.pkl"
    output:
        f"{RESULTS}/MALLET/merged_cistopic_with_models.pkl"
    params:
        mallet=config["paths"]["mallet_path"],
        n_topics=config["mallet"]["n_topics"],
        n_iter=config["mallet"]["n_iter"],
        n_cpu=config["resources"]["mallet_n_cpu"],
        tmp_path=config["mallet"]["tmp_path"],
        save_path=config["mallet"]["save_path"],
        alpha=config["mallet"]["alpha"],
        alpha_by_topic=config["mallet"]["alpha_by_topic"],
        eta=config["mallet"]["eta"],
        eta_by_topic=config["mallet"]["eta_by_topic"],
        random_state=config["mallet"]["random_state"],
        mallet_memory=config["mallet"]["mallet_memory"]
    shell:
        "python run_mallet_fixed.py "
        "--cistopic_obj_pickle {input} "
        "--mallet_path {params.mallet} "
        "--n_topics {params.n_topics} "
        "--n_cpu {params.n_cpu} "
        "--n_iter {params.n_iter} "
        "--tmp_path {params.tmp_path} "
        "--save_path {params.save_path} "
        "--mallet_memory {params.mallet_memory} "
        "--random_state {params.random_state} "
        "--alpha {params.alpha} "
        "{'--alpha_by_topic' if params.alpha_by_topic else ''} "
        "--eta {params.eta} "
        "{'--eta_by_topic' if params.eta_by_topic else ''}"

rule cluster_cistopic:
    input:
        f"{RESULTS}/MALLET/merged_cistopic_with_models.pkl"
    output:
        pickle=f"{RESULTS}/umap_clusters/cistopic_obj_clustered.pkl"
    params:
        ddir=f"{RESULTS}/umap_clusters",
        resolutions=" ".join(map(str, config["clustering"]["resolutions"])),
        k=config["clustering"]["k"]
    shell:
        "python cluster_cistopic_fixed.py "
        "-i {input} "
        "-o {output.pickle} "
        "-d {params.ddir} "
        "--resolutions {params.resolutions} "
        "--k {params.k}"

rule binarize_topics:
    input:
        pickle=f"{RESULTS}/umap_clusters/cistopic_obj_clustered.pkl"
    output:
        directory=f"{RESULTS}/topics"
    params:
        n=config["binarize"]["n"],
        s=config["binarize"]["s"]
    shell:
        "python binarize_topics_fixed.py "
        "--input_pickle {input.pickle} "
        "--output_dir {output.directory} "
        "-n {params.n} "
        "-s {params.s}"

rule dar_analysis:
    input:
        pickle=f"{RESULTS}/topics/cistopic_obj_binarized.pkl"
    output:
        directory=f"{RESULTS}/DAR_results"
    params:
        n_cpu=config["resources"]["dar_n_cpu"],
        scale_impute=config["dar"]["scale_impute"],
        scale_norm=config["dar"]["scale_norm"],
        adjpval_thr=config["dar"]["adjpval_thr"],
        log2fc_thr=config["dar"]["log2fc_thr"]
    shell:
        "python dar_analysis_fixed.py "
        "-i {input.pickle} "
        "-o {output.directory} "
        "-v celltype "
        "--n_cpu {params.n_cpu} "
        "--scale_impute {params.scale_impute} "
        "--scale_norm {params.scale_norm} "
        "--adjpval_thr {params.adjpval_thr} "
        "--log2fc_thr {params.log2fc_thr}"

rule export_region_sets:
    input:
        f"{RESULTS}/DAR_results/cistopic_obj_with_DARs.pkl"
    output:
        directory=RESULTS
    shell:
        "python export_region_sets_fixed.py -i {input} -o {output.directory}"

