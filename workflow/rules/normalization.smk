# Normalization and signal tracks: port of 1.3_cutrun_bigwig.sh,
# 1.7_cutrun_greenlist_sizefactors.sh and 1.8_cutrun_greenlist_bigwig.sh,
# plus spike-in scaling.
#
#   depth     : bamCoverage --normalizeUsing CPM | RPGC
#   greenlist : fragments at the CUT&RUN greenlist (de Mello et al. 2024;
#               multiBamSummary BED-file --extendReads --centerReads) ->
#               DESeq2 median-of-ratios size factors per normalization group
#   spikein   : spike-in fragments / their geometric mean per normalization group
#
# Size-factor tracks are scaled 1e6 / (size_factor x geometric-mean reads of
# the normalization group), i.e. CPM-equivalent: without a global shift they
# equal the CPM tracks, so all three methods share one unit. (Raw coverage /
# size factor would not be comparable with the CPM tracks.) Every track lives under results/bigwig/<method>/.

if "greenlist" in SIZEFACTOR_METHODS:

    rule greenlist_counts:
        input:
            bed=greenlist_bed(),
            bams=all_lib_bams(),
            bais=[f"{b}.bai" for b in all_lib_bams()],
        output:
            counts="results/normalization/greenlist/greenlist_counts.tsv",
        log:
            "logs/normalization/greenlist_counts.log",
        conda:
            "../envs/deeptools.yaml"
        threads: threads("deeptools", 8)
        resources:
            mem_mb=8000,
            runtime=240,
        params:
            labels=LIBS,
            mapq=config["align"]["min_mapq"],
            npz=lambda wildcards, output: output.counts.replace(".tsv", ".npz"),
        shell:
            """
            multiBamSummary BED-file --BED {input.bed} --bamfiles {input.bams} \
                --labels {params.labels} --extendReads --centerReads \
                --minMappingQuality {params.mapq} --numberOfProcessors {threads} \
                -o {params.npz} --outRawCounts {output.counts} > {log} 2>&1
            """


rule size_factors:
    input:
        flagstat=expand("results/qc/flagstat/{lib}.flagstat.txt", lib=LIBS),
        counts=lambda wildcards: (
            "results/normalization/greenlist/greenlist_counts.tsv"
            if wildcards.method == "greenlist"
            else "results/qc/spikein_summary.tsv"
        ),
    output:
        "results/normalization/{method}/size_factors.tsv",
    log:
        "logs/normalization/{method}.size_factors.log",
    wildcard_constraints:
        method="greenlist|spikein",
    conda:
        "../envs/python.yaml"
    params:
        method=lambda wildcards: wildcards.method,
        libs=LIBS,
        groups=[lib_group(l) for l in LIBS],
        targets=[LIBRARIES.loc[l, "target"] for l in LIBS],
        norm_groups=[norm_group(l) for l in LIBS],
    script:
        "../scripts/size_factors.py"


rule bigwig:
    input:
        bam="results/bam/{lib}.final.bam",
        bai="results/bam/{lib}.final.bam.bai",
        sf=lambda wildcards: (
            [] if wildcards.method == "depth" else sizefactor_table(wildcards.method)
        ),
        blacklist=blacklist_input(),
    output:
        "results/bigwig/{method}/{lib}.bw",
    log:
        "logs/bigwig/{method}/{lib}.log",
    conda:
        "../envs/deeptools.yaml"
    threads: threads("deeptools", 8)
    resources:
        mem_mb=8000,
        runtime=240,
    params:
        depth=NORM["depth_method"],
        gsize=REF["effective_genome_size"],
        bin=NORM["bin_size"],
        smooth=NORM["smooth_length"],
        mapq=config["align"]["min_mapq"],
        extra=NORM.get("bamcoverage_extra", ""),
    shell:
        """
        (
        set -euo pipefail
        if [ -z "{input.sf}" ]; then
            norm="--normalizeUsing {params.depth}"
            if [ "{params.depth}" = RPGC ]; then
                norm="$norm --effectiveGenomeSize {params.gsize}"
            fi
        else
            scale=$(awk -F'\\t' -v s={wildcards.lib} 'NR > 1 && $1 == s {{print $NF}}' {input.sf})
            [ -n "$scale" ] || {{ echo "no scale factor for {wildcards.lib}"; exit 1; }}
            norm="--scaleFactor $scale"
        fi
        bl=""
        if [ -n "{input.blacklist}" ]; then bl="--blackListFileName {input.blacklist}"; fi
        echo "bamCoverage $norm"
        bamCoverage --bam {input.bam} --outFileName {output} --outFileFormat bigwig \
            --extendReads --binSize {params.bin} --smoothLength {params.smooth} \
            --minMappingQuality {params.mapq} $norm $bl {params.extra} \
            --numberOfProcessors {threads}
        ) > {log} 2>&1
        """


rule group_bigwig:
    """Replicate-averaged track per group (bigwigAverage)."""
    input:
        lambda wildcards: [
            bigwig(wildcards.method, l) for l in LIBS_BY_GROUP[wildcards.group]
        ],
    output:
        "results/bigwig/{method}/groups/{group}.bw",
    log:
        "logs/bigwig/{method}/groups/{group}.log",
    conda:
        "../envs/deeptools.yaml"
    threads: threads("deeptools", 8)
    resources:
        mem_mb=8000,
        runtime=240,
    params:
        bin=NORM["bin_size"],
    shell:
        """
        (
        if [ $(echo {input} | wc -w) -eq 1 ]; then
            cp {input} {output}
        else
            bigwigAverage --bigwigs {input} --binSize {params.bin} \
                --numberOfProcessors {threads} --outFileName {output}
        fi
        ) > {log} 2>&1
        """
