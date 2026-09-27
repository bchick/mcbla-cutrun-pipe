# Peak calling: port of 1.4_cutrun_peaks.sh and 1.5_cutrun_consensus.sh,
# with IgG handling added.
#   1. merge replicate BAMs per group (and pool all IgG libraries)
#   2. MACS2 -f BAMPE --keep-dup all, per library and per merged group,
#      -q 0.05 (narrow) or --broad; with the group's IgG as -c when
#      igg.as_control is true
#   3. SEACR (opt-in) on fragment bedGraphs, IgG as control when available
#   4. blacklist filter -> <name>.raw.bed (narrowPeak columns)
#   5. IgG enrichment gate (igg.qc / igg.hotspot_filter) -> <name>.bed
#   6. per-group reproducible peaks (>= peaks.min_overlap replicates) and
#      per-target consensus / merged sets
# Downstream steps always read <name>.bed.


# ---------------------------------------------------------------------------
# Merged BAMs
# ---------------------------------------------------------------------------
rule merge_group_bam:
    input:
        merge_group_inputs,
    output:
        bam="results/bam/{kind}/{group}.bam",
        bai="results/bam/{kind}/{group}.bam.bai",
    log:
        "logs/peaks/{kind}/{group}.merge.log",
    wildcard_constraints:
        kind="merged|merged_sized",
    conda:
        "../envs/align.yaml"
    threads: threads("samtools", 8)
    resources:
        mem_mb=4000,
        runtime=240,
    shell:
        """
        (
        if [ $(echo {input} | wc -w) -eq 1 ]; then
            ln -sf "$(realpath {input})" {output.bam}
        else
            samtools merge -f -@ {threads} {output.bam} {input}
        fi
        samtools index -@ {threads} {output.bam}
        ) > {log} 2>&1
        """


rule igg_pool_bam:
    """All control libraries pooled (IgG hotspots; gate for groups without a control)."""
    input:
        igg_pool_inputs,
    output:
        bam="results/bam/igg_pool/pool{sfx}.bam",
        bai="results/bam/igg_pool/pool{sfx}.bam.bai",
    log:
        "logs/peaks/igg_pool{sfx}.merge.log",
    wildcard_constraints:
        sfx="|_sized",
    conda:
        "../envs/align.yaml"
    threads: threads("samtools", 8)
    shell:
        """
        (
        samtools merge -f -@ {threads} {output.bam} {input}
        samtools index -@ {threads} {output.bam}
        ) > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# MACS2
# ---------------------------------------------------------------------------
# narrowPeak (or broadPeak + a -1 summit column) -> blacklist -> raw.bed
MACS2_SHELL = """
(
set -euo pipefail
ctl=""
if [ -n "{input.control}" ]; then ctl="-c {input.control}"; fi
macs2 callpeak -t {input.bam} $ctl -f BAMPE -g {params.gsize} --keep-dup all \
    {params.mode} --outdir {params.outdir} --name {params.name} {params.extra}
native={params.outdir}/{params.name}_peaks.narrowPeak
if [ -e {params.outdir}/{params.name}_peaks.broadPeak ]; then
    native={params.outdir}/{params.name}_peaks.broadPeak
fi
awk 'BEGIN{{OFS="\\t"}} {{ if (NF < 10) $10 = -1; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10 }}' \
    "$native" > {output.peaks}.tmp
if [ -n "{input.blacklist}" ]; then
    bedtools intersect -v -a {output.peaks}.tmp -b {input.blacklist} > {output.peaks}
else
    mv {output.peaks}.tmp {output.peaks}
fi
rm -f {output.peaks}.tmp
echo "peaks: $(wc -l < "$native") called, $(wc -l < {output.peaks}) after blacklist"
) > {log} 2>&1
"""


rule macs2:
    input:
        bam=lambda wildcards: peak_bam(wildcards.name),
        control=lambda wildcards: control_bam(wildcards.name),
        blacklist=blacklist_input(),
    output:
        peaks="results/peaks/macs2/{level}/{name}.raw.bed",
        xls="results/peaks/macs2/{level}/{name}_peaks.xls",
    log:
        "logs/peaks/macs2/{level}/{name}.log",
    conda:
        "../envs/macs2.yaml"
    resources:
        mem_mb=lambda wildcards: 8000 if wildcards.level == "individual" else 16000,
        runtime=480,
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.peaks),
        name=lambda wildcards: wildcards.name,
        mode=lambda wildcards: macs2_mode_args(wildcards.name),
        gsize=REF["macs2_gsize"],
        extra=PK.get("macs2_extra", ""),
    shell:
        MACS2_SHELL


rule macs2_relaxed:
    """Relaxed (-p 0.01) narrow peaks: IDR input (libraries) and oracle (groups)."""
    input:
        bam=lambda wildcards: peak_bam(wildcards.name),
        control=lambda wildcards: control_bam(wildcards.name),
        blacklist=blacklist_input(),
    output:
        peaks="results/peaks/idr/relaxed/{name}.raw.bed",
        xls="results/peaks/idr/relaxed/{name}_peaks.xls",
    log:
        "logs/peaks/idr/relaxed/{name}.log",
    conda:
        "../envs/macs2.yaml"
    resources:
        mem_mb=16000,
        runtime=480,
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.peaks),
        name=lambda wildcards: wildcards.name,
        mode=lambda wildcards: macs2_mode_args(wildcards.name, relaxed=True),
        gsize=REF["macs2_gsize"],
        extra=PK.get("macs2_extra", ""),
    shell:
        MACS2_SHELL


# ---------------------------------------------------------------------------
# SEACR (opt-in)
# ---------------------------------------------------------------------------
if PK["seacr"]["run"]:

    rule fragment_bedgraph:
        """Fragment coverage (< 1 kb, same-contig pairs) for SEACR."""
        input:
            bam="{prefix}.bam",
            sizes="results/reference/host.chrom.sizes",
        output:
            temp("{prefix}.fragments.bedgraph"),
        log:
            "logs/seacr/bedgraph/{prefix}.log",
        wildcard_constraints:
            prefix="results/bam/.+",
        conda:
            "../envs/align.yaml"
        threads: threads("samtools", 4)
        resources:
            mem_mb=16000,
            runtime=240,
        params:
            tmp=lambda wildcards: f"{wildcards.prefix}.bg_tmp",
        shell:
            """
            (
            set -euo pipefail
            samtools sort -n -@ {threads} -T {params.tmp} -o {params.tmp}.n.bam {input.bam}
            bedtools bamtobed -bedpe -i {params.tmp}.n.bam 2> /dev/null \
              | awk 'BEGIN{{OFS="\\t"}} $1 == $4 && $6 - $2 < 1000 {{print $1, $2, $6}}' \
              | LC_ALL=C sort -k1,1 -k2,2n -S 4G -T $(dirname {params.tmp}) \
              | bedtools genomecov -bg -i - -g {input.sizes} > {output}
            rm -f {params.tmp}.n.bam
            ) > {log} 2>&1
            """

    rule seacr:
        input:
            bg=lambda wildcards: peak_bam(wildcards.name).replace(
                ".bam", ".fragments.bedgraph"
            ),
            control=lambda wildcards: [
                c.replace(".bam", ".fragments.bedgraph")
                for c in control_bam(wildcards.name)
            ],
            blacklist=blacklist_input(),
        output:
            peaks="results/peaks/seacr/{level}/{name}.raw.bed",
        log:
            "logs/peaks/seacr/{level}/{name}.log",
        conda:
            "../envs/seacr.yaml"
        resources:
            mem_mb=16000,
            runtime=240,
        params:
            prefix=lambda wildcards, output: output.peaks.replace(".raw.bed", ""),
            threshold=PK["seacr"]["threshold"],
            norm=PK["seacr"]["norm"],
            mode=PK["seacr"]["mode"],
            name=lambda wildcards: wildcards.name,
        shell:
            """
            (
            set -euo pipefail
            ctl="{params.threshold}"
            if [ -n "{input.control}" ]; then ctl="{input.control}"; fi
            SEACR_1.3.sh {input.bg} $ctl {params.norm} {params.mode} {params.prefix}
            # chr start end total max max_region -> narrowPeak (summit = max region centre)
            awk -v n={params.name} 'BEGIN{{OFS="\\t"}} {{
                    split($6, r, /[:-]/); s = int((r[2] + r[3]) / 2) - $2
                    sc = int($4); if (sc > 1000) sc = 1000
                    print $1, $2, $3, n "_seacr_" NR, sc, ".", $4, -1, -1, s }}' \
                {params.prefix}.{params.mode}.bed > {output.peaks}.tmp
            if [ -n "{input.blacklist}" ]; then
                bedtools intersect -v -a {output.peaks}.tmp -b {input.blacklist} > {output.peaks}
                rm -f {output.peaks}.tmp
            else
                mv {output.peaks}.tmp {output.peaks}
            fi
            echo "SEACR peaks after blacklist: $(wc -l < {output.peaks})"
            ) > {log} 2>&1
            """


# ---------------------------------------------------------------------------
# IgG hotspots, background windows and the enrichment gate
# ---------------------------------------------------------------------------
if USE_HOTSPOTS:

    rule igg_hotspots:
        """Pooled-IgG MACS2 peaks (plus igg.hotspot_bed), merged."""
        input:
            bam="results/bam/igg_pool/pool.bam" if HAS_IGG else [],
            extra=[IGG["hotspot_bed"]] if IGG.get("hotspot_bed") else [],
        output:
            "results/igg/hotspots.bed",
        log:
            "logs/igg/hotspots.log",
        conda:
            "../envs/macs2.yaml"
        resources:
            mem_mb=16000,
            runtime=240,
        params:
            outdir="results/igg/macs2",
            q=IGG["hotspot_qvalue"],
            gsize=REF["macs2_gsize"],
        shell:
            """
            (
            set -euo pipefail
            : > {output}.tmp
            if [ -n "{input.bam}" ]; then
                macs2 callpeak -t {input.bam} -f BAMPE -g {params.gsize} --keep-dup all \
                    -q {params.q} --outdir {params.outdir} --name igg_pool
                cut -f1-3 {params.outdir}/igg_pool_peaks.narrowPeak >> {output}.tmp
            fi
            for b in {input.extra}; do cut -f1-3 "$b" >> {output}.tmp; done
            sort -k1,1 -k2,2n {output}.tmp | bedtools merge > {output}
            rm -f {output}.tmp
            echo "hotspots: $(wc -l < {output})"
            ) > {log} 2>&1
            """


if IGG_FOLD:

    rule background_windows:
        """Random windows away from peaks and the blacklist (per-library depth scale).

        Enrichment is depth-free: each library is divided by its own mean
        density over these windows (the design of the mcf7 analysis 39 IgG
        reality gate), not by total reads, which a high-FRiP target inflates.
        """
        input:
            sizes="results/reference/host.chrom.sizes",
            peaks=expand(
                "results/peaks/{caller}/merged/{g}.raw.bed",
                caller=CALLERS,
                g=TARGET_GROUPS,
            ),
            blacklist=blacklist_input(),
        output:
            "results/igg/background_windows.bed",
        log:
            "logs/igg/background_windows.log",
        conda:
            "../envs/align.yaml"
        params:
            n=IGG["background_windows"],
            w=IGG["window_size"],
            seed=IGG["seed"],
        shell:
            """
            (
            set -euo pipefail
            tmp=$(mktemp -d)
            awk -v w={params.w} '$2 >= 10 * w' {input.sizes} > $tmp/genome
            cat {input.peaks} {input.blacklist} | cut -f1-3 \
              | awk 'BEGIN{{OFS="\\t"}} {{s = $2 - 1000; if (s < 0) s = 0; print $1, s, $3 + 1000}}' \
              | sort -k1,1 -k2,2n | bedtools merge > $tmp/excl.bed
            bedtools random -l {params.w} -n {params.n} -seed {params.seed} -g $tmp/genome \
              | bedtools shuffle -i - -g $tmp/genome -excl $tmp/excl.bed -noOverlapping \
                  -seed {params.seed} -maxTries 10000 \
              | cut -f1-3 | sort -k1,1 -k2,2n > {output}
            rm -rf $tmp
            echo "background windows: $(wc -l < {output})"
            ) > {log} 2>&1
            """

    rule igg_counts:
        """Fragment centres of target and IgG in each peak and background window."""
        input:
            peaks="results/peaks/{caller}/{level}/{name}.raw.bed",
            bg="results/igg/background_windows.bed",
            target=lambda wildcards: peak_bam(wildcards.name),
            igg=lambda wildcards: igg_reference_bam(wildcards.name),
            hotspots="results/igg/hotspots.bed",
        output:
            counts=temp("results/peaks/{caller}/{level}/{name}.igg_counts.tsv"),
            hot=temp("results/peaks/{caller}/{level}/{name}.hotspot_overlap.tsv"),
        log:
            "logs/igg/{caller}/{level}/{name}.counts.log",
        conda:
            "../envs/deeptools.yaml"
        threads: threads("deeptools", 8)
        resources:
            mem_mb=8000,
            runtime=240,
        params:
            mapq=config["align"]["min_mapq"],
            tmp=lambda wildcards, output: f"{output.counts}.work",
        shell:
            """
            (
            set -euo pipefail
            mkdir -p {params.tmp}
            if [ -s {input.peaks} ]; then
                cut -f1-3 {input.peaks} > {params.tmp}/peaks.bed
                cut -f1-3 {input.bg} > {params.tmp}/bg.bed
                multiBamSummary BED-file --BED {params.tmp}/peaks.bed {params.tmp}/bg.bed \
                    --bamfiles {input.target} {input.igg} --labels target igg \
                    --extendReads --centerReads --minMappingQuality {params.mapq} \
                    --numberOfProcessors {threads} \
                    -o {params.tmp}/counts.npz --outRawCounts {output.counts}
            else
                printf "#'chr'\\t'start'\\t'end'\\t'target'\\t'igg'\\n" > {output.counts}
            fi
            cut -f1-3 {input.peaks} | bedtools intersect -c -a - -b {input.hotspots} \
                > {output.hot} 2> /dev/null || : > {output.hot}
            rm -rf {params.tmp}
            ) > {log} 2>&1
            """

    rule igg_gate:
        """Per-peak target / IgG enrichment; filter when igg.hotspot_filter is on."""
        input:
            peaks="results/peaks/{caller}/{level}/{name}.raw.bed",
            bg="results/igg/background_windows.bed",
            counts="results/peaks/{caller}/{level}/{name}.igg_counts.tsv",
            hot="results/peaks/{caller}/{level}/{name}.hotspot_overlap.tsv",
        output:
            peaks="results/peaks/{caller}/{level}/{name}.bed",
            table="results/peaks/{caller}/{level}/{name}.igg.tsv",
            summary="results/peaks/{caller}/{level}/{name}.igg_summary.tsv",
        log:
            "logs/igg/{caller}/{level}/{name}.gate.log",
        conda:
            "../envs/python.yaml"
        params:
            name=lambda wildcards: wildcards.name,
            min_fold=IGG["min_fold"],
            filter=IGG_FILTER,
            remove_hotspots=bool(IGG["remove_hotspot_overlaps"]) and IGG_FILTER,
        script:
            "../scripts/igg_gate.py"

else:

    rule finalize_peaks:
        """No IgG gate: raw peaks, minus igg.hotspot_bed overlaps when filtering."""
        input:
            peaks="results/peaks/{caller}/{level}/{name}.raw.bed",
            hotspots="results/igg/hotspots.bed" if IGG_FILTER else [],
        output:
            peaks="results/peaks/{caller}/{level}/{name}.bed",
        log:
            "logs/peaks/{caller}/{level}/{name}.finalize.log",
        conda:
            "../envs/align.yaml"
        shell:
            """
            (
            if [ -n "{input.hotspots}" ]; then
                bedtools intersect -v -a {input.peaks} -b {input.hotspots} > {output.peaks}
            else
                cp {input.peaks} {output.peaks}
            fi
            ) > {log} 2>&1
            """


rule igg_enrichment_table:
    """Per-library and per-group IgG enrichment summaries (MultiQC, qc_summary)."""
    input:
        ind=expand(
            "results/peaks/{caller}/individual/{lib}.igg_summary.tsv",
            caller=PRIMARY,
            lib=TARGET_LIBS,
        ),
        grp=expand(
            "results/peaks/{caller}/merged/{g}.igg_summary.tsv",
            caller=PRIMARY,
            g=TARGET_GROUPS,
        ),
    output:
        "results/qc/igg_enrichment.tsv",
    log:
        "logs/igg/igg_enrichment_table.log",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (head -n 1 {input.ind[0]}; for f in {input.ind} {input.grp}; do tail -n +2 $f; done) \
            > {output} 2> {log}
        """


# ---------------------------------------------------------------------------
# Peak statistics
# ---------------------------------------------------------------------------
rule peak_stats:
    """Peaks, FRiP (fragments in peaks / fragments, duplicates excluded), median width."""
    input:
        peaks="results/peaks/{caller}/{level}/{name}.bed",
        bam=lambda wildcards: (
            lib_bam(wildcards.name)
            if wildcards.level == "individual"
            else f"results/bam/merged/{wildcards.name}.bam"
        ),
    output:
        temp("results/peaks/{caller}/{level}/{name}.stats.tsv"),
    log:
        "logs/peaks/{caller}/{level}/{name}.stats.log",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (
        set -euo pipefail
        n=$(wc -l < {input.peaks})
        frip=NA; median=NA
        if [ "$n" -gt 0 ]; then
            total=$(samtools view -c -f 66 -F 3332 {input.bam})
            inpk=$(samtools view -c -f 66 -F 3332 -L {input.peaks} {input.bam})
            frip=$(awk -v a="$inpk" -v b="$total" 'BEGIN{{ if (b>0) printf "%.4f", a/b; else print "NA" }}')
            median=$(awk '{{print $3-$2}}' {input.peaks} | sort -n \
              | awk '{{a[n++]=$1}} END{{ if (n%2==1) print a[int(n/2)]; else print (a[n/2-1]+a[n/2])/2 }}')
        fi
        printf "%s\\t%s\\t%s\\t%s\\n" "{wildcards.name}" "$n" "$frip" "$median" > {output}
        ) > {log} 2>&1
        """


rule peak_summary:
    input:
        lambda wildcards: expand(
            "results/peaks/{caller}/{level}/{name}.stats.tsv",
            caller=wildcards.caller,
            level=wildcards.level,
            name=peak_names(wildcards.level),
        ),
    output:
        "results/peaks/{caller}/{level}/peak_summary.tsv",
    log:
        "logs/peaks/{caller}/{level}/peak_summary.log",
    conda:
        "../envs/align.yaml"
    params:
        label=lambda wildcards: (
            "Sample" if wildcards.level == "individual" else "Group"
        ),
    shell:
        """
        (printf "{params.label}\\tPeaks\\tFRiP\\tMedian_Width\\n"; cat {input}) > {output} 2> {log}
        """


rule frip_library:
    """FRiP of one library over its group's merged-BAM peaks (primary caller)."""
    input:
        bam="results/bam/{lib}.final.bam",
        bai="results/bam/{lib}.final.bam.bai",
        peaks=lambda wildcards: frip_peaks(wildcards.lib),
    output:
        "results/qc/frip/{lib}.frip.tsv",
    log:
        "logs/qc/{lib}.frip.log",
    conda:
        "../envs/align.yaml"
    threads: threads("samtools", 4)
    shell:
        """
        (
        set -euo pipefail
        total=$(samtools view -@ {threads} -c -f 66 -F 3332 {input.bam})
        inpk=0
        if [ -s {input.peaks} ]; then
            inpk=$(samtools view -@ {threads} -c -f 66 -F 3332 -L {input.peaks} {input.bam})
        fi
        frip=$(awk -v a="$inpk" -v b="$total" 'BEGIN{{ if (b > 0) printf "%.4f", a / b; else print "NA" }}')
        printf "Sample\\tFragments\\tFragments_In_Peaks\\tFRiP\\n%s\\t%s\\t%s\\t%s\\n" \
            "{wildcards.lib}" "$total" "$inpk" "$frip" > {output}
        ) > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Consensus peak sets
# ---------------------------------------------------------------------------
rule group_consensus:
    """Peaks present in >= min(peaks.min_overlap, n replicates) libraries of a group.

    Replicate peaks are pooled and merged; a merged region is kept when peaks
    from at least k distinct replicates fall in it (full peak extents are kept,
    unlike an intersection). Grouping is by samplesheet group, so arms such as
    vehicle / inhibitor at the same time point are never pooled.
    """
    input:
        lambda wildcards: [
            peaks_final(PRIMARY, "individual", l)
            for l in LIBS_BY_GROUP[wildcards.group]
        ],
    output:
        "results/peaks/consensus/groups/{group}.bed",
    log:
        "logs/peaks/consensus/{group}.log",
    conda:
        "../envs/align.yaml"
    params:
        k=lambda wildcards: group_rep_threshold(wildcards.group),
    shell:
        """
        (
        set -euo pipefail
        i=0
        for f in {input}; do
            i=$((i + 1)); awk -v i=$i 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, i}}' "$f"
        done | sort -k1,1 -k2,2n > {output}.tmp
        if [ -s {output}.tmp ]; then
            bedtools merge -i {output}.tmp -c 4 -o count_distinct \
              | awk -v k={params.k} 'BEGIN{{OFS="\\t"}} $4 >= k {{print $1, $2, $3}}' > {output}
        else
            : > {output}
        fi
        rm -f {output}.tmp
        echo "{wildcards.group}: $(wc -l < {output}) peaks in >= {params.k} replicates"
        ) > {log} 2>&1
        """


rule target_peakset:
    """Union over a target's groups (bedtools merge)."""
    input:
        target_set_inputs,
    output:
        "results/peaks/consensus/{target}.{kind}.bed",
    log:
        "logs/peaks/consensus/{target}.{kind}.log",
    wildcard_constraints:
        kind="consensus|merged|idr",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (cat {input} | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > {output}) 2> {log}
        """


rule consensus_summary:
    input:
        groups=expand("results/peaks/consensus/groups/{g}.bed", g=TARGET_GROUPS),
        targets=expand(
            "results/peaks/consensus/{t}.{k}.bed", t=TARGETS, k=["consensus", "merged"]
        ),
    output:
        "results/peaks/consensus/consensus_summary.tsv",
    log:
        "logs/peaks/consensus/summary.log",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (
        printf "Set\\tPeaks\\tMedian_Width\\n"
        for f in {input.groups} {input.targets}; do
            n=$(wc -l < $f)
            med=$(awk '{{print $3-$2}}' $f | sort -n \
              | awk '{{a[n++]=$1}} END{{ if (n==0) print "NA"; else if (n%2==1) print a[int(n/2)]; else print (a[n/2-1]+a[n/2])/2 }}')
            printf "%s\\t%s\\t%s\\n" "$(basename $f .bed)" "$n" "$med"
        done
        ) > {output} 2> {log}
        """
