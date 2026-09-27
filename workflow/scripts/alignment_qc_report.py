"""Assemble the per-library alignment QC table (port of the QC TSV row written
by 1.1_cutrun_align_cc.sh / 1.1_cutrun_align_spikein.sh).

Columns: Sample, Raw_Reads, Trimmed_Reads, Aligned_Reads, Aligned_Pct,
Mito_Reads, Mito_Pct, Spikein_Reads_Raw, Blacklist_Removed, Final_Reads,
Dup_Pct, Mean_FragSize, TF_Fraction, NRF, PBC1, PBC2

Differences from the shell original (deliberate):
  * Raw_Reads comes from cutadapt's "Total read pairs processed" instead of a
    separate zcat pass over R1 (same number, one fewer full read of the FASTQ).
  * Dup_Pct is DUPLICATE TOTAL / EXAMINED from `samtools markdup -f`.
  * Final_Reads is the flagstat total (= `samtools view -c` on the final BAM;
    marked duplicates included for targets, removed for controls).
  * Fragment sizes exclude duplicates; TF_Fraction is the fraction of
    fragments <= fragments.tf_max_size (default 120 bp).
"""

import re
import sys

sm = snakemake  # noqa: F821  (injected by Snakemake)
sys.stderr = open(sm.log[0], "w")


def grab(path, pattern, cast=int, default="NA"):
    with open(path) as fh:
        text = fh.read()
    m = re.search(pattern, text, flags=re.MULTILINE)
    if not m:
        return default
    return cast(m.group(1).replace(",", ""))


def kv(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                out[parts[0]] = parts[1]
    return out


def pct(num, den):
    try:
        return f"{100.0 * float(num) / float(den):.1f}%"
    except (ValueError, ZeroDivisionError, TypeError):
        return "NA"


header = [
    "Sample", "Raw_Reads", "Trimmed_Reads", "Aligned_Reads", "Aligned_Pct",
    "Mito_Reads", "Mito_Pct", "Spikein_Reads_Raw", "Blacklist_Removed",
    "Final_Reads", "Dup_Pct", "Mean_FragSize", "TF_Fraction", "NRF", "PBC1", "PBC2",
]
tf_max = int(sm.params.tf_max)

rows = []
for i, lib in enumerate(sm.params.libs):
    cut = sm.input.cutadapt[i]
    bt2 = sm.input.bowtie2[i]
    raw = grab(cut, r"Total read pairs processed:\s+([\d,]+)")
    trimmed = grab(cut, r"Pairs written \(passing filters\):\s+([\d,]+)")
    aligned_pct = grab(bt2, r"([\d.]+)% overall alignment rate", cast=str)
    filt = kv(sm.input.filt[i])
    total = filt.get("total_aligned", "NA")
    mito = filt.get("mito_reads", "NA")
    spike = filt.get("spikein_reads_raw", "NA")
    bl = filt.get("blacklist_removed", "NA")
    final = grab(sm.input.flagstat[i], r"^(\d+) \+ \d+ in total")
    examined = grab(sm.input.markdup[i], r"^EXAMINED:\s+(\d+)")
    dups = grab(sm.input.markdup[i], r"^DUPLICATE TOTAL:\s+(\d+)")
    s = n = small = 0
    with open(sm.input.frag[i]) as fh:
        for line in fh:
            size, count = (int(x) for x in line.split())
            s += size * count
            n += count
            if size <= tf_max:
                small += count
    mean_frag = f"{s / n:.0f}" if n else "NA"
    tf_frac = f"{small / n:.3f}" if n else "NA"
    with open(sm.input.complexity[i]) as fh:
        cx = dict(zip(*(line.rstrip("\n").split("\t") for line in fh)))
    rows.append([
        lib, raw, trimmed, total,
        f"{aligned_pct}%" if aligned_pct != "NA" else "NA",
        mito, pct(mito, total), spike, bl, final, pct(dups, examined), mean_frag,
        tf_frac, cx.get("NRF", "NA"), cx.get("PBC1", "NA"), cx.get("PBC2", "NA"),
    ])

with open(sm.output[0], "w") as out:
    out.write("\t".join(header) + "\n")
    for r in rows:
        out.write("\t".join(str(x) for x in r) + "\n")
