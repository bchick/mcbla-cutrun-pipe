# IDR (opt-in, modules.idr): the ATAC module (mcbla-bulkatac-pipe idr.smk,
# a port of 2.2_atac_idr.sh) applied to CUT&RUN narrow groups.
#
#   relaxed MACS2 (-p 0.01, with the IgG control when used) per replicate and
#   per merged group (oracle) -> sort -k8,8rn -> idr --rank p.value
#   --peak-list <oracle> --idr-threshold 0.05 -> keep column 5 >= 540
#
# Replicates: 2 = one pair; > 2 = every pair, keeping the pair with the most
# reproducible peaks (encode_max_pair) or rep1 vs rep2 (first_two).
# Single-replicate and broad groups use their merged peaks in the idr set.

if MODULES.get("idr", False):

    rule sort_relaxed_peaks:
        input:
            "results/peaks/idr/relaxed/{name}.raw.bed",
        output:
            temp("results/peaks/idr/sorted/{name}.sorted.narrowPeak"),
        log:
            "logs/idr/{name}.sort.log",
        conda:
            "../envs/align.yaml"
        shell:
            "sort -k8,8rn {input} > {output} 2> {log}"

    rule idr_pair:
        input:
            a="results/peaks/idr/sorted/{a}.sorted.narrowPeak",
            b="results/peaks/idr/sorted/{b}.sorted.narrowPeak",
            oracle="results/peaks/idr/sorted/{group}.sorted.narrowPeak",
        output:
            all="results/peaks/idr/pairs/{group}/{a}__{b}.idr_all.narrowPeak",
            filt="results/peaks/idr/pairs/{group}/{a}__{b}.idr.narrowPeak",
            status="results/peaks/idr/pairs/{group}/{a}__{b}.status",
        log:
            "logs/idr/{group}/{a}__{b}.idr.log",
        wildcard_constraints:
            a=_alt(LIBS),
            b=_alt(LIBS),
        conda:
            "../envs/idr.yaml"
        resources:
            mem_mb=8000,
            runtime=120,
        params:
            thr=IDR["threshold"],
            scaled=idr_scaled_threshold(),
            rank=IDR.get("rank", "p.value"),
            allow_failure="1" if IDR.get("allow_failure", True) else "0",
        shell:
            """
            if idr --samples {input.a} {input.b} --peak-list {input.oracle} \
                   --input-file-type narrowPeak --rank {params.rank} \
                   --output-file {output.all} --idr-threshold {params.thr} \
                   --plot > {log} 2>&1; then
                awk -v t={params.scaled} 'BEGIN{{OFS="\\t"}} $5 >= t' {output.all} \
                    | cut -f1-10 > {output.filt}
                echo OK > {output.status}
            elif [ "{params.allow_failure}" = 1 ]; then
                echo "WARNING: IDR failed for {wildcards.a} vs {wildcards.b}; recorded as FAILED" >> {log}
                : > {output.all}; : > {output.filt}
                echo FAILED > {output.status}
            else
                exit 1
            fi
            """

    rule idr_select:
        """Pick the reproducible set per group (max pair when > 2 replicates)."""
        input:
            pairs=idr_pair_files,
            oracle="results/peaks/idr/relaxed/{group}.raw.bed",
            reps=lambda wildcards: expand(
                "results/peaks/idr/relaxed/{lib}.raw.bed",
                lib=LIBS_BY_GROUP[wildcards.group],
            ),
        output:
            peaks="results/peaks/idr/{group}_idr.narrowPeak",
            stats=temp("results/peaks/idr/{group}.idr_stats.tsv"),
        log:
            "logs/idr/{group}.select.log",
        conda:
            "../envs/python.yaml"
        params:
            condition=lambda wildcards: wildcards.group,
            libs=lambda wildcards: LIBS_BY_GROUP[wildcards.group],
            strategy=IDR.get("multi_rep_strategy", "encode_max_pair"),
        script:
            "../scripts/idr_select.py"

    rule idr_summary:
        input:
            stats=expand("results/peaks/idr/{g}.idr_stats.tsv", g=IDR_GROUPS),
        output:
            "results/peaks/idr/idr_summary.tsv",
        log:
            "logs/idr/idr_summary.log",
        conda:
            "../envs/python.yaml"
        params:
            single=[],
            fallback="stringent",
        script:
            "../scripts/idr_summary.py"
