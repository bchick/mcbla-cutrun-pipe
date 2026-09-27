"""Spike-in fragments per library (fastq mode: counted after alignment; bam
mode: the samplesheet `spikein_reads` column) and as a percentage of host
fragments (mapped reads / 2 in the final BAM)."""

import csv
import re
import sys

sm = snakemake  # noqa: F821
sys.stderr = open(sm.log[0], "w")

libs = list(sm.params.libs)
spike = {}
for path in sm.input.get("counts", []):
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            spike[row["Sample"]] = int(row["Spikein_Fragments"])
if not spike:
    spike = {lib: int(v) for lib, v in zip(libs, sm.params.from_sheet) if v != ""}

with open(sm.output[0], "w") as out:
    out.write("Sample\tHost_Fragments\tSpikein_Fragments\tSpikein_Pct\n")
    for lib, path in zip(libs, sm.input.flagstat):
        with open(path) as fh:
            m = re.search(r"^(\d+) \+ \d+ mapped", fh.read(), flags=re.MULTILINE)
        host = int(m.group(1)) // 2 if m else 0
        s = spike.get(lib)
        pct = f"{100 * s / host:.3f}" if s is not None and host else "NA"
        out.write(f"{lib}\t{host}\t{s if s is not None else 'NA'}\t{pct}\n")
