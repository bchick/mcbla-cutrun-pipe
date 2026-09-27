# QC: port of 1.2_cutrun_qc.sh (deepTools) plus FastQC, a per-library QC
# summary with CUT&RUN thresholds, and MultiQC.
#   bamPEFragmentSize, multiBamSummary 500 bp bins (blacklist excluded) ->
#   Spearman heatmap + PCA, plotFingerprint (with the pooled IgG as the
#   JS-distance reference when there are IgG libraries).

DT = config["qc"]

if INPUT_MODE == "fastq":

    rule fastqc:
        """FastQC on the raw reads (runs of one library concatenated)."""
        input:
            unpack(trim_inputs),
        output:
            r1="results/qc/fastqc/{lib}_R1_fastqc.zip",
            r2="results/qc/fastqc/{lib}_R2_fastqc.zip",
        log:
            "logs/qc/{lib}.fastqc.log",
        conda:
            "../envs/fastqc.yaml"
        threads: 2
        resources:
            mem_mb=2000,
            runtime=240,
        shell:
            """
            (
            set -euo pipefail
            tmp=$(mktemp -d)
            # name the reads after the library so MultiQC shows library names
            ln -s "$(realpath {input.r1})" "$tmp/{wildcards.lib}_R1.fastq.gz"
            ln -s "$(realpath {input.r2})" "$tmp/{wildcards.lib}_R2.fastq.gz"
            fastqc --threads {threads} --outdir "$tmp" \
                "$tmp/{wildcards.lib}_R1.fastq.gz" "$tmp/{wildcards.lib}_R2.fastq.gz"
            mv "$tmp/{wildcards.lib}_R1_fastqc.zip" {output.r1}
            mv "$tmp/{wildcards.lib}_R2_fastqc.zip" {output.r2}
            rm -rf "$tmp"
            ) > {log} 2>&1
            """


rule fragment_size_plot:
    input:
        bam=all_lib_bams(),
        bai=[f"{b}.bai" for b in all_lib_bams()],
    output:
        png="results/qc/deeptools/fragment_size_distribution.png",
        table="results/qc/deeptools/fragment_size_table.tsv",
        raw="results/qc/deeptools/fragment_size_raw.tsv",
    log:
        "logs/qc/bamPEFragmentSize.log",
    conda:
        "../envs/deeptools.yaml"
    threads: threads("deeptools", 8)
    resources:
        mem_mb=8000,
        runtime=240,
    params:
        labels=LIBS,
    shell:
        """
        bamPEFragmentSize --bamfiles {input.bam} --samplesLabel {params.labels} \
            --numberOfProcessors {threads} --maxFragmentLength 1000 \
            --plotTitle "Fragment Size Distribution (CUT&RUN)" \
            --table {output.table} --outRawFragmentLengths {output.raw} \
            -o {output.png} > {log} 2>&1
        """


rule multibamsummary:
    input:
        bam=all_lib_bams(),
        bai=[f"{b}.bai" for b in all_lib_bams()],
        blacklist=blacklist_input(),
    output:
        "results/qc/deeptools/multiBamSummary.npz",
    log:
        "logs/qc/multiBamSummary.log",
    conda:
        "../envs/deeptools.yaml"
    threads: threads("deeptools", 8)
    resources:
        mem_mb=16000,
        runtime=480,
    params:
        labels=LIBS,
        bin=DT["summary_bin_size"],
        mapq=config["align"]["min_mapq"],
        bl=lambda wildcards, input: (
            f"--blackListFileName {input.blacklist}" if input.blacklist else ""
        ),
    shell:
        """
        multiBamSummary bins --bamfiles {input.bam} --labels {params.labels} \
            --binSize {params.bin} --extendReads --minMappingQuality {params.mapq} \
            {params.bl} --numberOfProcessors {threads} \
            --outFileName {output} > {log} 2>&1
        """


rule plot_correlation:
    input:
        "results/qc/deeptools/multiBamSummary.npz",
    output:
        png="results/qc/deeptools/correlation_spearman.png",
        tsv="results/qc/deeptools/correlation_spearman.tsv",
    log:
        "logs/qc/plotCorrelation.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        plotCorrelation --corData {input} --corMethod spearman --whatToPlot heatmap \
            --plotNumbers --skipZeros --plotTitle "Spearman Correlation (CUT&RUN)" \
            --outFileCorMatrix {output.tsv} --plotFile {output.png} > {log} 2>&1
        """


rule plot_pca:
    input:
        "results/qc/deeptools/multiBamSummary.npz",
    output:
        png="results/qc/deeptools/pca_plot.png",
        tsv="results/qc/deeptools/pca_data.tsv",
    log:
        "logs/qc/plotPCA.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        plotPCA --corData {input} --plotTitle "PCA of CUT&RUN libraries" \
            --outFileNameData {output.tsv} --plotFile {output.png} > {log} 2>&1
        """


rule fingerprint:
    input:
        bam=all_lib_bams(),
        bai=[f"{b}.bai" for b in all_lib_bams()],
        igg="results/bam/igg_pool/pool.bam" if HAS_IGG else [],
    output:
        png="results/qc/deeptools/fingerprint.png",
        metrics="results/qc/deeptools/fingerprint_metrics.tsv",
        counts="results/qc/deeptools/fingerprint_counts.tsv",
    log:
        "logs/qc/plotFingerprint.log",
    conda:
        "../envs/deeptools.yaml"
    threads: threads("deeptools", 8)
    resources:
        mem_mb=8000,
        runtime=240,
    params:
        # the pooled IgG is plotted too, as the JS-distance reference
        labels=LIBS + (["IgG_pool"] if HAS_IGG else []),
        jsd=lambda wildcards, input: f"--JSDsample {input.igg}" if input.igg else "",
    shell:
        """
        plotFingerprint --bamfiles {input.bam} {input.igg} --labels {params.labels} \
            --extendReads --numberOfProcessors {threads} {params.jsd} \
            --plotTitle "Fingerprint (CUT&RUN)" \
            --outQualityMetrics {output.metrics} --outRawCounts {output.counts} \
            --plotFile {output.png} > {log} 2>&1
        """


rule qc_summary:
    input:
        unpack(qc_summary_inputs),
    output:
        "results/qc/qc_summary.tsv",
    log:
        "logs/qc/qc_summary.log",
    conda:
        "../envs/python.yaml"
    params:
        libs=LIBS,
        groups=[lib_group(l) for l in LIBS],
        roles=["control" if is_control_lib(l) else "target" for l in LIBS],
        types=[LIBRARIES.loc[l, "target_type"] for l in LIBS],
        tf_max=config["fragments"]["tf_max_size"],
        thresholds=QC_THRESHOLDS,
    script:
        "../scripts/qc_summary.py"


rule multiqc:
    input:
        multiqc_inputs,
    output:
        "results/qc/multiqc/multiqc_report.html",
    log:
        "logs/qc/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        config=os.path.join(workflow.basedir, "resources", "multiqc_config.yaml"),
    shell:
        """
        (
        stage=$(mktemp -d)
        for f in {input}; do
            base=$(basename "$f")
            case "$f" in
              *alignment_qc_report.tsv) dest=alignment_qc_mqc.tsv ;;
              *macs2/individual/peak_summary.tsv) dest=peaks_macs2_individual_mqc.tsv ;;
              *macs2/merged/peak_summary.tsv) dest=peaks_macs2_merged_mqc.tsv ;;
              *seacr/individual/peak_summary.tsv) dest=peaks_seacr_individual_mqc.tsv ;;
              *seacr/merged/peak_summary.tsv) dest=peaks_seacr_merged_mqc.tsv ;;
              *idr_summary.tsv) dest=idr_summary_mqc.tsv ;;
              *qc_summary.tsv) dest=qc_summary_mqc.tsv ;;
              *igg_enrichment.tsv) dest=igg_enrichment_mqc.tsv ;;
              *spikein_summary.tsv) dest=spikein_mqc.tsv ;;
              *greenlist/size_factors.tsv) dest=sizefactors_greenlist_mqc.tsv ;;
              *spikein/size_factors.tsv) dest=sizefactors_spikein_mqc.tsv ;;
              *) dest="$base" ;;
            esac
            ln -sf "$(realpath "$f")" "$stage/$dest"
        done
        multiqc --force --config {params.config} -o {params.outdir} "$stage"
        rm -rf "$stage"
        ) > {log} 2>&1
        """
