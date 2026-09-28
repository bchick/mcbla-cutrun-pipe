# Differential binding (opt-in, diff.run): DiffBind with the DESeq2 backend,
# one DiffBind object per target (contrasts compare groups of one target).
#
#   sample sheet -> dba.count over the diff.peaks set (or minOverlap
#   consensus for `individual`) -> dba.normalize under diff.normalization:
#     depth     DiffBind default (library size)
#     greenlist dba.normalize(library = greenlist size factors, DBA_NORM_LIB)
#     spikein   dba.normalize(library = spike-in size factors, DBA_NORM_LIB)
#     csaw      background bins (normcheck only)
#   -> dba.contrast (~group, or ~batch + group) -> dba.analyze -> dba.report
#
# The greenlist / spike-in route follows the lab's greenlist helper for
# DiffBind. normcheck (opt-in) re-runs the same contrasts on
# the same counted object under the other methods and flags contrasts whose
# gained/lost counts move by more than normcheck.threshold.

if RUN_DIFF:

    rule diffbind_samplesheet:
        input:
            bams=lambda wildcards: all_lib_bams(diff_libs(wildcards.target)),
            peaks=lambda wildcards: [
                peaks_final(PRIMARY, "individual", l)
                for l in diff_libs(wildcards.target)
            ],
        output:
            "results/diff/{target}/diffbind_samplesheet.csv",
        log:
            "logs/diff/{target}/diffbind_samplesheet.log",
        conda:
            "../envs/python.yaml"
        params:
            libs=lambda wildcards: diff_libs(wildcards.target),
            groups=lambda wildcards: [lib_group(l) for l in diff_libs(wildcards.target)],
            replicates=lambda wildcards: [
                LIBRARIES.loc[l, "replicate"] for l in diff_libs(wildcards.target)
            ],
            # DiffBind "Factor" carries the batch when config `batch` is true
            factors=lambda wildcards: [
                LIBRARIES.loc[l, "batch"] if USE_BATCH else wildcards.target
                for l in diff_libs(wildcards.target)
            ],
            tissue=lambda wildcards: wildcards.target,
        script:
            "../scripts/diffbind_samplesheet.py"

    rule diffbind_count:
        input:
            sheet="results/diff/{target}/diffbind_samplesheet.csv",
            bams=lambda wildcards: all_lib_bams(diff_libs(wildcards.target)),
            consensus=lambda wildcards: DIFF_PEAKS[1](wildcards.target) or [],
        output:
            "results/diff/{target}/dba_counted.rds",
        log:
            "logs/diff/{target}/diffbind_count.log",
        conda:
            "../envs/r.yaml"
        threads: threads("diffbind", 8)
        resources:
            mem_mb=64000,
            runtime=720,
        params:
            min_overlap=DIFF["min_overlap"],
            summits=DIFF.get("summits"),
            peakset=DIFF_PEAKS[0],
        script:
            "../scripts/diffbind_count.R"

    rule diffbind_analyze:
        input:
            dba="results/diff/{target}/dba_counted.rds",
            contrasts=config["contrasts"],
            sf=lambda wildcards: (
                sizefactor_table(wildcards.dmethod)
                if wildcards.dmethod in ("greenlist", "spikein")
                else []
            ),
        output:
            rds="results/diff/{target}/{dmethod}/dba_analyzed.rds",
            summary="results/diff/{target}/{dmethod}/summary.tsv",
            sizefactors="results/diff/{target}/{dmethod}/size_factors.tsv",
            tables=directory("results/diff/{target}/{dmethod}/tables"),
            plots="results/diff/{target}/{dmethod}/plots.pdf",
        log:
            "logs/diff/{target}/{dmethod}.analyze.log",
        wildcard_constraints:
            dmethod="depth|greenlist|spikein|csaw",
        conda:
            "../envs/r.yaml"
        threads: threads("diffbind", 8)
        resources:
            mem_mb=32000,
            runtime=480,
        params:
            helpers=os.path.join(workflow.basedir, "scripts", "common.R"),
            method=lambda wildcards: wildcards.dmethod,
            labels=lambda wildcards: target_contrasts(wildcards.target),
            batch=USE_BATCH,
            fdr=DIFF["fdr"],
            lfc=DIFF["lfc"],
        script:
            "../scripts/diffbind_analyze.R"

    rule contrast_beds:
        """Gained / lost peaks of one contrast as BED (heatmaps, annotation, motifs)."""
        input:
            tables=lambda wildcards: (
                f"results/diff/{CONTRAST_TARGET[wildcards.label]}/{DIFF_NORM}/tables"
            ),
        output:
            gained="results/diff/contrasts/{label}_gained.bed",
            lost="results/diff/contrasts/{label}_lost.bed",
        log:
            "logs/diff/contrasts/{label}.beds.log",
        conda:
            "../envs/python.yaml"
        params:
            fdr=DIFF["fdr"],
            lfc=DIFF["lfc"],
        script:
            "../scripts/contrast_beds.py"


if RUN_NORMCHECK:

    rule normcheck_compare:
        input:
            primary="results/diff/{target}/%s/summary.tsv" % DIFF_NORM,
            others=expand("results/diff/{{target}}/{m}/summary.tsv", m=NC_METHODS),
        output:
            table="results/diff/{target}/normcheck/norm_comparison.tsv",
            verdict="results/diff/{target}/normcheck/norm_verdict.tsv",
            barplot="results/diff/{target}/normcheck/norm_comparison_barplot.pdf",
        log:
            "logs/diff/{target}/normcheck.log",
        conda:
            "../envs/r.yaml"
        params:
            helpers=os.path.join(workflow.basedir, "scripts", "common.R"),
            primary=DIFF_NORM,
            methods=NC_METHODS,
            labels=lambda wildcards: target_contrasts(wildcards.target),
            fdr=DIFF["fdr"],
            threshold=NC["sensitivity_threshold"],
            min_abs=NC["min_abs_change"],
        script:
            "../scripts/normcheck_compare.R"
