# Downstream analyses (opt-in), each over the peak set named in <section>.peaks.
#
#   heatmaps : computeMatrix reference-point --referencePoint center +/- flank
#              on replicate-averaged bigWigs of each target's groups (and their
#              IgG controls) over the target's peaks, TSSs and, with diff on,
#              each contrast's gained / lost peaks -> plotHeatmap
#              (settings of the lab's heatmap scripts)
#   annotate : ChIPseeker annotatePeak with a TxDb built from reference.gtf
#              (genome-agnostic; gene names from the GTF)
#   motifs   : HOMER findMotifsGenome.pl on reference.fasta, known motifs
#              (plus de novo with motifs.denovo); target peaks against a
#              GC-matched genomic background, contrast peaks against the
#              target's peak set

HM = config["heatmaps"]
ANN = config["annotate"]
MOT = config["motifs"]


# Placeholder figure when a region set is empty (e.g. no significant peaks).
EMPTY_PNG = """python - <<'PY'
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(4, 2))
fig.text(0.5, 0.5, "{msg}", ha="center", va="center")
fig.savefig("{out}")
PY"""


rule tss_sites:
    """Unique transcript TSSs from the GTF (gene TSSs if it has no transcripts)."""
    input:
        gtf=REF["gtf"] if REF.get("gtf") else [],
    output:
        "results/reference/tss.bed",
    log:
        "logs/reference/tss_sites.log",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (
        set -euo pipefail
        for feature in transcript gene; do
            zcat -f {input.gtf} \
              | awk -F'\\t' -v f=$feature 'BEGIN{{OFS="\\t"}} $3 == f {{
                    p = ($7 == "-") ? $5 - 1 : $4 - 1; print $1, p, p + 1, ".", ".", $7}}' \
              | LC_ALL=C sort -u -k1,1 -k2,2n -k6,6 > {output}
            [ -s {output} ] && break
        done
        echo "TSSs: $(wc -l < {output})"
        ) > {log} 2>&1
        """


if RUN_HEATMAPS:
    HM_NORM = HM.get("normalization", "depth")

    rule heatmap_target:
        input:
            regions=heatmap_regions,
            bw=lambda wildcards: [
                group_bigwig(HM_NORM, g) for g in heatmap_groups(wildcards.target)
            ],
        output:
            matrix="results/heatmaps/{target}/{kind}_matrix.gz",
            png="results/heatmaps/{target}/{kind}_heatmap.png",
        log:
            "logs/heatmaps/{target}.{kind}.log",
        wildcard_constraints:
            kind="peaks|tss",
        conda:
            "../envs/deeptools.yaml"
        threads: threads("deeptools", 8)
        resources:
            mem_mb=16000,
            runtime=240,
        params:
            labels=lambda wildcards: heatmap_groups(wildcards.target),
            flank=HM["flank"],
            ref=lambda wildcards: "center" if wildcards.kind == "peaks" else "TSS",
            cmap=HM.get("colormap", "Blues"),
            title=lambda wildcards: f"{wildcards.target} ({HM_NORM} normalization)",
            empty=lambda wildcards, output: EMPTY_PNG.format(
                msg="no regions", out=output.png
            ),
        shell:
            """
            (
            set -euo pipefail
            if [ ! -s {input.regions} ]; then
                : | gzip > {output.matrix}
                {params.empty}
                exit 0
            fi
            computeMatrix reference-point --referencePoint {params.ref} \
                -S {input.bw} --samplesLabel {params.labels} -R {input.regions} \
                -a {params.flank} -b {params.flank} --missingDataAsZero --skipZeros \
                -p {threads} -o {output.matrix}
            plotHeatmap -m {output.matrix} -o {output.png} --colorMap {params.cmap} \
                --heatmapHeight 20 --heatmapWidth 3 --plotTitle "{params.title}"
            ) > {log} 2>&1
            """

    if RUN_DIFF:

        rule heatmap_contrast:
            input:
                gained="results/diff/contrasts/{label}_gained.bed",
                lost="results/diff/contrasts/{label}_lost.bed",
                bw=lambda wildcards: [
                    group_bigwig(HM_NORM, g) for g in CONTRAST_GROUPS[wildcards.label]
                ],
            output:
                png="results/heatmaps/contrasts/{label}_heatmap.png",
            log:
                "logs/heatmaps/contrasts/{label}.log",
            conda:
                "../envs/deeptools.yaml"
            threads: threads("deeptools", 8)
            resources:
                mem_mb=16000,
                runtime=240,
            params:
                labels=lambda wildcards: list(CONTRAST_GROUPS[wildcards.label]),
                flank=HM["flank"],
                cmap=HM.get("colormap", "Blues"),
                matrix=lambda wildcards, output: output.png.replace(
                    ".png", "_matrix.gz"
                ),
                empty=lambda wildcards, output: EMPTY_PNG.format(
                    msg="no differential peaks", out=output.png
                ),
            shell:
                """
                (
                set -euo pipefail
                regions=""; names=""
                for f in {input.gained} {input.lost}; do
                    if [ -s "$f" ]; then
                        regions="$regions $f"
                        names="$names $(basename $f .bed | sed 's/^{wildcards.label}_//')"
                    fi
                done
                if [ -z "$regions" ]; then
                    {params.empty}
                    exit 0
                fi
                computeMatrix reference-point --referencePoint center \
                    -S {input.bw} --samplesLabel {params.labels} -R $regions \
                    -a {params.flank} -b {params.flank} --missingDataAsZero --skipZeros \
                    -p {threads} -o {params.matrix}
                plotHeatmap -m {params.matrix} -o {output.png} --colorMap {params.cmap} \
                    --regionsLabel $names --heatmapHeight 20 --heatmapWidth 3 \
                    --plotTitle "{wildcards.label}"
                ) > {log} 2>&1
                """


if RUN_ANNOTATE:
    ANN_SETS = {t: ANN_PEAKS[1](t) for t in TARGETS}
    ANN_SETS.update(contrast_sets())

    rule annotate_peaks:
        input:
            beds=list(ANN_SETS.values()),
            gtf=REF["gtf"],
        output:
            summary="results/annotate/annotation_summary.tsv",
            plots="results/annotate/annotation_plots.pdf",
            tables=expand("results/annotate/{name}.annotation.tsv", name=list(ANN_SETS)),
        log:
            "logs/annotate/annotate.log",
        conda:
            "../envs/annotate.yaml"
        resources:
            mem_mb=16000,
            runtime=240,
        params:
            names=list(ANN_SETS),
            tss=ANN["tss_region"],
        script:
            "../scripts/annotate_peaks.R"


# HOMER on one foreground BED; -bg only for contrast sets. A private
# preparsed directory per job avoids concurrent jobs racing on one cache.
HOMER_SHELL = """
(
set -euo pipefail
mkdir -p {params.outdir}
if [ ! -s {input.bed} ]; then
    echo "no peaks in {input.bed}" > {output}
    exit 0
fi
work=$(mktemp -d -p {params.outdir})
to_homer() {{ awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, "peak" NR, ".", "+"}}' "$1"; }}
to_homer {input.bed} > $work/fg.bed
bg=""
if [ -n "{params.bg}" ]; then
    to_homer {params.bg} > $work/bg.bed
    bg="-bg $work/bg.bed"
fi
findMotifsGenome.pl $work/fg.bed {input.fasta} {params.outdir} \\
    -size {params.size} -p {threads} {params.denovo} {params.mask} $bg \\
    -preparsedDir $work/preparsed {params.extra}
rm -rf $work
# HOMER skips the known-motif step (exit 0, no output) for very small sets
if [ ! -s {output} ]; then
    echo "HOMER wrote no known-motif results ($(wc -l < {input.bed}) peaks; too few)" > {output}
fi
) > {log} 2>&1
"""

if RUN_MOTIFS:
    HOMER_PARAMS = dict(
        size=MOT["size"],
        denovo="" if MOT.get("denovo") else "-nomotif",
        mask="-mask" if MOT.get("mask") else "",
        extra=MOT.get("extra", ""),
    )

    rule homer_target:
        input:
            bed=lambda wildcards: MOTIF_PEAKS[1](wildcards.target),
            fasta="results/reference/genome.fa",
        output:
            "results/motifs/{target}/knownResults.txt",
        log:
            "logs/motifs/{target}.log",
        conda:
            "../envs/homer.yaml"
        threads: threads("homer", 8)
        resources:
            mem_mb=16000,
            runtime=720,
        params:
            **HOMER_PARAMS,
            outdir=lambda wildcards, output: os.path.dirname(output[0]),
            bg="",
        shell:
            HOMER_SHELL

    rule homer_contrast:
        input:
            bed="results/diff/contrasts/{label}_{direction}.bed",
            bg=lambda wildcards: MOTIF_PEAKS[1](CONTRAST_TARGET[wildcards.label]),
            fasta="results/reference/genome.fa",
        output:
            "results/motifs/contrasts/{label}_{direction}/knownResults.txt",
        log:
            "logs/motifs/contrasts/{label}_{direction}.log",
        wildcard_constraints:
            direction="gained|lost",
        conda:
            "../envs/homer.yaml"
        threads: threads("homer", 8)
        resources:
            mem_mb=16000,
            runtime=720,
        params:
            **HOMER_PARAMS,
            outdir=lambda wildcards, output: os.path.dirname(output[0]),
            bg=lambda wildcards, input: input.bg,
        shell:
            HOMER_SHELL
