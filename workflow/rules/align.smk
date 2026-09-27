# Alignment and filtering (input_mode: fastq), or import of filtered BAMs
# (input_mode: bam).
#
# FASTQ mode is a port of the lab's 1.1_cutrun_align_cc.sh and
# 1.1_cutrun_align_spikein.sh (mcf7_project/scripts/alignment/):
#   cutadapt (NEBNext/TruSeq, -q 20,20 -m 20 --pair-filter=any)
#   -> bowtie2 --very-sensitive-local --no-mixed --no-discordant -I 10 -X 700 --dovetail
#      (host, or host + spike-in combined index)
#   -> samtools view -q 20 -f 2 -F 2828; spike-in reads split off and counted
#   -> drop mito, bedtools intersect -v blacklist (optional)
#   -> sort -n -> fixmate -m -> sort -> markdup -s
#      (targets: duplicates marked and kept, the lab convention;
#       controls: removed, as nf-core/cutandrun does)
#   -> index, flagstat, fragment sizes, alignment QC table.


# ---------------------------------------------------------------------------
# Reference preparation
# ---------------------------------------------------------------------------
rule prepare_fasta:
    input:
        fasta=REF["fasta"],
    output:
        fasta="results/reference/genome.fa",
        fai="results/reference/genome.fa.fai",
        sizes="results/reference/chrom.sizes",
    log:
        "logs/reference/prepare_fasta.log",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (
        case "{input.fasta}" in
            *.gz) pigz -dc "{input.fasta}" > {output.fasta} ;;
            *)    ln -sf "$(realpath {input.fasta})" {output.fasta} ;;
        esac
        samtools faidx {output.fasta}
        cut -f1,2 {output.fai} > {output.sizes}
        ) > {log} 2>&1
        """


rule host_chrom_sizes:
    """Host contigs without the mito contig (background windows, shuffles)."""
    input:
        "results/reference/chrom.sizes",
    output:
        "results/reference/host.chrom.sizes",
    log:
        "logs/reference/host_chrom_sizes.log",
    conda:
        "../envs/align.yaml"
    params:
        mito=MITO,
        prefix=SPIKE_PREFIX,
    shell:
        """
        awk -v m={params.mito} -v p={params.prefix} \
            '$1 != m && index($1, p) != 1' {input} > {output} 2> {log}
        """


if INPUT_MODE == "fastq" and SPIKE_ALIGN and not SPIKE.get("combined_bowtie2_index"):

    rule combined_fasta:
        """Host + spike-in FASTA; spike-in contigs renamed <contig_prefix><name>."""
        input:
            host="results/reference/genome.fa",
            spike=SPIKE["fasta"],
        output:
            "results/reference/combined.fa",
        log:
            "logs/reference/combined_fasta.log",
        conda:
            "../envs/align.yaml"
        params:
            prefix=SPIKE_PREFIX,
        shell:
            """
            (
            set -euo pipefail
            cat {input.host} > {output}
            zcat -f {input.spike} | awk -v p={params.prefix} \
                'substr($0, 1, 1) == ">" {{print ">" p substr($1, 2); next}} {{print}}' >> {output}
            grep -c "^>{params.prefix}" {output}
            ) > {log} 2>&1
            """


if INPUT_MODE == "fastq":

    rule bowtie2_build:
        input:
            fasta=lambda wildcards: (
                "results/reference/combined.fa"
                if wildcards.which == "bowtie2_combined"
                else "results/reference/genome.fa"
            ),
        output:
            multiext(
                "results/reference/{which}/genome",
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ),
        log:
            "logs/reference/{which}_build.log",
        wildcard_constraints:
            which="bowtie2|bowtie2_combined",
        conda:
            "../envs/align.yaml"
        threads: threads("bowtie2_build", 8)
        resources:
            mem_mb=16000,
            runtime=480,
        params:
            prefix=lambda wildcards, output: output[0][: -len(".1.bt2")],
        shell:
            "bowtie2-build --threads {threads} {input.fasta} {params.prefix} > {log} 2>&1"

    rule merge_fastq:
        """Concatenate sequencing runs of one library (nf-core semantics)."""
        input:
            merge_fastq_inputs,
        output:
            temp("results/fastq/{lib}_R{read}.merged.fastq.gz"),
        log:
            "logs/align/{lib}_R{read}.merge_fastq.log",
        wildcard_constraints:
            read="1|2",
        conda:
            "../envs/align.yaml"
        shell:
            "cat {input} > {output} 2> {log}"

    rule cutadapt:
        input:
            unpack(trim_inputs),
        output:
            r1=temp("results/fastq/{lib}_R1.trimmed.fastq.gz"),
            r2=temp("results/fastq/{lib}_R2.trimmed.fastq.gz"),
        log:
            "logs/align/{lib}.cutadapt.log",
        conda:
            "../envs/align.yaml"
        threads: threads("cutadapt", 8)
        resources:
            mem_mb=4000,
            runtime=240,
        params:
            a1=config["trimming"]["adapter_r1"],
            a2=config["trimming"]["adapter_r2"],
            q=config["trimming"]["quality"],
            m=config["trimming"]["min_length"],
            extra=config["trimming"].get("extra", ""),
        shell:
            """
            cutadapt -a {params.a1} -A {params.a2} \
                -q {params.q} -m {params.m} --pair-filter=any \
                -j {threads} {params.extra} \
                -o {output.r1} -p {output.r2} \
                {input.r1} {input.r2} > {log} 2>&1
            """

    rule bowtie2_align:
        input:
            r1="results/fastq/{lib}_R1.trimmed.fastq.gz",
            r2="results/fastq/{lib}_R2.trimmed.fastq.gz",
            idx=bowtie2_index_files(),
        output:
            bam=temp("results/bam/tmp/{lib}.aligned.bam"),
        log:
            "logs/align/{lib}.bowtie2.log",
        conda:
            "../envs/align.yaml"
        threads: threads("bowtie2", 16)
        resources:
            mem_mb=16000,
            runtime=720,
        params:
            index=lambda wildcards, input: re.sub(r"\.1\.bt2l?$", "", input.idx[0]),
            args=config["align"]["bowtie2_args"],
        shell:
            """
            (set -o pipefail
            bowtie2 -p {threads} -x {params.index} {params.args} \
                -1 {input.r1} -2 {input.r2} 2> {log} \
            | samtools view -b -@ 2 -o {output.bam} -) 2>> {log}
            """

    rule filter_bam:
        """MAPQ / proper-pair / flag filter; split off spike-in; drop mito; blacklist."""
        input:
            bam="results/bam/tmp/{lib}.aligned.bam",
            blacklist=blacklist_input(),
        output:
            bam=temp("results/bam/tmp/{lib}.filtered.bam"),
            spike=temp("results/bam/tmp/{lib}.spikein.bam") if SPIKE_ALIGN else [],
            stats="results/qc/filter_stats/{lib}.filter_stats.tsv",
        log:
            "logs/align/{lib}.filter_bam.log",
        conda:
            "../envs/align.yaml"
        threads: threads("samtools", 8)
        resources:
            mem_mb=8000,
            runtime=240,
        params:
            mapq=config["align"]["min_mapq"],
            flag_exclude=config["align"]["flag_exclude"],
            mito=MITO,
            spike_bam=lambda wildcards, output: output.spike if SPIKE_ALIGN else "",
            prefix=SPIKE_PREFIX,
        shell:
            """
            (
            set -euo pipefail
            tmp={output.bam}.host.bam
            # total, mito and spike-in records in the raw alignment
            samtools view -@ {threads} {input.bam} \
              | awk -v m={params.mito} -v p={params.prefix} \
                  'BEGIN{{n=0;c=0;s=0}} {{n++}} $3==m{{c++}} index($3,p)==1{{s++}}
                   END{{print n"\\t"c"\\t"s}}' > {output.stats}.counts
            read total mito spike < {output.stats}.counts
            rm -f {output.stats}.counts

            samtools view -h -@ {threads} -q {params.mapq} -f 2 -F {params.flag_exclude} \
                -o {output.bam}.q.bam {input.bam}
            if [ -n "{params.spike_bam}" ]; then
                # both mates on spike-in contigs
                samtools view -h {output.bam}.q.bam \
                  | awk -v p={params.prefix} '$1 ~ /^@/ || (index($3,p)==1 && ($7=="=" || index($7,p)==1))' \
                  | samtools view -b -o {params.spike_bam} -
            fi
            # host stream: drop mito and spike-in contigs (either mate)
            samtools view -h {output.bam}.q.bam \
              | awk -v m={params.mito} -v p={params.prefix} \
                  '$1 ~ /^@/ || ($3 != m && index($3,p)!=1 && index($7,p)!=1)' \
              | samtools view -b -@ 2 -o $tmp -
            rm -f {output.bam}.q.bam
            post_mito=$(samtools view -c $tmp)

            if [ -n "{input.blacklist}" ]; then
                bedtools intersect -v -abam $tmp -b {input.blacklist} > {output.bam}
                rm -f $tmp
            else
                mv $tmp {output.bam}
            fi
            post_bl=$(samtools view -c {output.bam})

            printf "total_aligned\\t%s\\nmito_reads\\t%s\\nspikein_reads_raw\\t%s\\npost_mito_filter\\t%s\\npost_blacklist\\t%s\\nblacklist_removed\\t%s\\n" \
                "$total" "$mito" "$spike" "$post_mito" "$post_bl" "$((post_mito - post_bl))" > {output.stats}
            ) > {log} 2>&1
            """

    rule library_complexity:
        """ENCODE NRF / PBC1 / PBC2 on the filtered BAM before duplicate marking.

        One fragment per pair (the leftmost mate, TLEN > 0), keyed by
        chrom, start, TLEN and strand:
          NRF  = distinct fragments / all fragments
          PBC1 = fragments seen exactly once / distinct fragments
          PBC2 = fragments seen exactly once / fragments seen exactly twice
        """
        input:
            "results/bam/tmp/{lib}.filtered.bam",
        output:
            "results/qc/complexity/{lib}.complexity.tsv",
        log:
            "logs/align/{lib}.complexity.log",
        conda:
            "../envs/align.yaml"
        threads: threads("samtools", 4)
        resources:
            mem_mb=8000,
            runtime=240,
        params:
            tmp=lambda wildcards: f"results/bam/tmp/{wildcards.lib}.complexity",
        shell:
            """
            (
            set -euo pipefail
            mkdir -p {params.tmp}
            samtools view -@ {threads} -F 2308 {input} \
              | awk 'BEGIN{{OFS="\\t"}} $9 > 0 {{print $3, $4, $9, int($2 / 16) % 2}}' \
              | LC_ALL=C sort -S 4G --parallel={threads} -T {params.tmp} \
              | uniq -c \
              | awk 'BEGIN{{OFS="\\t"}}
                     {{t += $1; d++; if ($1 == 1) m1++; if ($1 == 2) m2++}}
                     END{{
                       nrf  = (t > 0) ? sprintf("%.4f", d / t) : "NA"
                       pbc1 = (d > 0) ? sprintf("%.4f", m1 / d) : "NA"
                       pbc2 = (m2 > 0) ? sprintf("%.4f", m1 / m2) : "NA"
                       print "Total_Fragments", "Distinct_Fragments", "One_Read", "Two_Reads", "NRF", "PBC1", "PBC2"
                       print t + 0, d + 0, m1 + 0, m2 + 0, nrf, pbc1, pbc2
                     }}' > {output}
            rm -rf {params.tmp}
            ) > {log} 2>&1
            """

    rule markdup:
        """sort -n -> fixmate -m -> sort -> markdup -s (-r for controls by default)."""
        input:
            "results/bam/tmp/{lib}.filtered.bam",
        output:
            bam="results/bam/{lib}.final.bam",
            stats="results/qc/markdup/{lib}.markdup.txt",
        log:
            "logs/align/{lib}.markdup.log",
        conda:
            "../envs/align.yaml"
        threads: threads("samtools", 8)
        resources:
            mem_mb=16000,
            runtime=240,
        params:
            tmp=lambda wildcards: f"results/bam/tmp/{wildcards.lib}",
            mem=config["align"].get("sort_mem_per_thread", "768M"),
            dedup=lambda wildcards: "-r" if remove_dups(wildcards.lib) else "",
        shell:
            """
            (
            set -euo pipefail
            samtools sort -@ {threads} -m {params.mem} -n -T {params.tmp}.nsort \
                -o {params.tmp}.nsorted.bam {input}
            samtools fixmate -@ {threads} -m {params.tmp}.nsorted.bam {params.tmp}.fixmate.bam
            rm -f {params.tmp}.nsorted.bam
            samtools sort -@ {threads} -m {params.mem} -T {params.tmp}.csort \
                -o {params.tmp}.csorted.bam {params.tmp}.fixmate.bam
            rm -f {params.tmp}.fixmate.bam
            samtools markdup -@ {threads} {params.dedup} -s -f {output.stats} \
                {params.tmp}.csorted.bam {output.bam}
            rm -f {params.tmp}.csorted.bam
            ) > {log} 2>&1
            """

    rule spikein_counts:
        """Deduplicated spike-in fragments (read 1 of proper pairs)."""
        input:
            "results/bam/tmp/{lib}.spikein.bam",
        output:
            "results/qc/spikein/{lib}.spikein.tsv",
        log:
            "logs/align/{lib}.spikein.log",
        conda:
            "../envs/align.yaml"
        threads: threads("samtools", 4)
        params:
            tmp=lambda wildcards: f"results/bam/tmp/{wildcards.lib}.spk",
        shell:
            """
            (
            set -euo pipefail
            raw=$(samtools view -c -f 64 {input})
            samtools sort -n -@ {threads} -T {params.tmp}.n -o {params.tmp}.n.bam {input}
            samtools fixmate -m {params.tmp}.n.bam {params.tmp}.fm.bam
            samtools sort -@ {threads} -T {params.tmp}.c -o {params.tmp}.c.bam {params.tmp}.fm.bam
            samtools markdup -r -s {params.tmp}.c.bam {params.tmp}.dd.bam
            dedup=$(samtools view -c -f 64 {params.tmp}.dd.bam)
            rm -f {params.tmp}.n.bam {params.tmp}.fm.bam {params.tmp}.c.bam {params.tmp}.dd.bam
            printf "Sample\\tSpikein_Fragments_Raw\\tSpikein_Fragments\\n%s\\t%s\\t%s\\n" \
                "{wildcards.lib}" "$raw" "$dedup" > {output}
            ) > {log} 2>&1
            """

    rule index_bam:
        input:
            "results/bam/{lib}.final.bam",
        output:
            "results/bam/{lib}.final.bam.bai",
        log:
            "logs/align/{lib}.index.log",
        conda:
            "../envs/align.yaml"
        threads: threads("samtools", 4)
        shell:
            "samtools index -@ {threads} {input} > {log} 2>&1"

    rule alignment_qc_report:
        input:
            cutadapt=expand("logs/align/{lib}.cutadapt.log", lib=LIBS),
            bowtie2=expand("logs/align/{lib}.bowtie2.log", lib=LIBS),
            filt=expand("results/qc/filter_stats/{lib}.filter_stats.tsv", lib=LIBS),
            markdup=expand("results/qc/markdup/{lib}.markdup.txt", lib=LIBS),
            complexity=expand("results/qc/complexity/{lib}.complexity.tsv", lib=LIBS),
            flagstat=expand("results/qc/flagstat/{lib}.flagstat.txt", lib=LIBS),
            frag=expand("results/qc/fragment_sizes/{lib}_fragment_sizes.tsv", lib=LIBS),
        output:
            "results/qc/alignment_qc_report.tsv",
        log:
            "logs/align/alignment_qc_report.log",
        conda:
            "../envs/python.yaml"
        params:
            libs=LIBS,
            tf_max=config["fragments"]["tf_max_size"],
        script:
            "../scripts/alignment_qc_report.py"


# ---------------------------------------------------------------------------
# BAM mode: filtered BAMs from an earlier run (e.g. the mcf7 *_markdup.bam)
# ---------------------------------------------------------------------------
if INPUT_MODE == "bam":

    rule import_bam:
        """Link the BAM; drop marked duplicates where duplicates.remove_* asks."""
        input:
            bam=import_bam_input,
        output:
            bam="results/bam/{lib}.final.bam",
            bai="results/bam/{lib}.final.bam.bai",
        log:
            "logs/import/{lib}.import.log",
        conda:
            "../envs/align.yaml"
        threads: threads("samtools", 4)
        params:
            dedup=lambda wildcards: "1" if remove_dups(wildcards.lib) else "0",
        shell:
            """
            (
            set -euo pipefail
            if [ "{params.dedup}" = 1 ]; then
                samtools view -b -@ {threads} -F 1024 -o {output.bam} {input.bam}
                samtools index -@ {threads} {output.bam}
            else
                ln -sf "$(realpath {input.bam})" {output.bam}
                if [ -s "{input.bam}.bai" ]; then
                    ln -sf "$(realpath {input.bam}.bai)" {output.bai}
                else
                    samtools index -@ {threads} {output.bam}
                fi
            fi
            ) > {log} 2>&1
            """


# ---------------------------------------------------------------------------
# Both modes
# ---------------------------------------------------------------------------
rule flagstat:
    input:
        bam="results/bam/{lib}.final.bam",
        bai="results/bam/{lib}.final.bam.bai",
    output:
        "results/qc/flagstat/{lib}.flagstat.txt",
    log:
        "logs/align/{lib}.flagstat.log",
    conda:
        "../envs/align.yaml"
    threads: threads("samtools", 4)
    shell:
        "samtools flagstat -@ {threads} {input.bam} > {output} 2> {log}"


rule fragment_sizes:
    """Fragment-length histogram (TLEN of read 1, duplicates excluded)."""
    input:
        bam="results/bam/{lib}.final.bam",
        bai="results/bam/{lib}.final.bam.bai",
    output:
        "results/qc/fragment_sizes/{lib}_fragment_sizes.tsv",
    log:
        "logs/align/{lib}.fragment_sizes.log",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (set -o pipefail
        samtools view -f 66 -F 3332 {input.bam} \
          | awk '$9 > 0 && $9 < 1000 {{print $9}}' \
          | sort -n | uniq -c | awk '{{print $2"\\t"$1}}' > {output}) 2> {log}
        """


rule size_select:
    """Fragments <= fragments.peak_max_size (TF peak calling on sub-nucleosomal fragments)."""
    input:
        bam="results/bam/{lib}.final.bam",
        bai="results/bam/{lib}.final.bam.bai",
    output:
        bam="results/bam/sized/{lib}.bam",
        bai="results/bam/sized/{lib}.bam.bai",
    log:
        "logs/align/{lib}.size_select.log",
    conda:
        "../envs/align.yaml"
    threads: threads("samtools", 4)
    params:
        max=PEAK_MAX_SIZE,
    shell:
        """
        (
        set -euo pipefail
        samtools view -h -@ {threads} {input.bam} \
          | awk -v m={params.max} '$1 ~ /^@/ || ($9 != 0 && $9 <= m && $9 >= -m)' \
          | samtools view -b -@ 2 -o {output.bam} -
        samtools index {output.bam}
        echo "kept $(samtools view -c {output.bam}) of $(samtools view -c {input.bam}) records"
        ) > {log} 2>&1
        """


rule spikein_summary:
    """Spike-in fragments per library, and as a fraction of host fragments."""
    input:
        flagstat=expand("results/qc/flagstat/{lib}.flagstat.txt", lib=LIBS),
        counts=spikein_count_files(),
    output:
        "results/qc/spikein_summary.tsv",
    log:
        "logs/qc/spikein_summary.log",
    conda:
        "../envs/python.yaml"
    params:
        libs=LIBS,
        from_sheet=[LIBRARIES.loc[l, "spikein_reads"] for l in LIBS],
    script:
        "../scripts/spikein_summary.py"
