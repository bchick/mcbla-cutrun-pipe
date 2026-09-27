"""Gained and lost peaks of one contrast (FDR < diff.fdr) as BED6.

Name = peak_id, score = round(-10 log10 FDR) capped at 1000, strand ".",
sorted by position. Empty files when nothing is significant.
"""

import csv
import math
import os
import sys

sm = snakemake  # noqa: F821
sys.stderr = open(sm.log[0], "w")

label = sm.wildcards.label
fdr = float(sm.params.fdr)
table = os.path.join(sm.input.tables, f"{label}_all.tsv")
gained, lost = [], []
with open(table) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        try:
            q, fold = float(row["FDR"]), float(row["Fold"])
        except ValueError:
            continue
        if q >= fdr:
            continue
        score = min(1000, int(round(-10 * math.log10(max(q, 1e-300)))))
        rec = (row["seqnames"], int(row["start"]) - 1, int(row["end"]), row["peak_id"], score)
        (gained if fold > 0 else lost).append(rec)

for recs, path in ((gained, sm.output.gained), (lost, sm.output.lost)):
    with open(path, "w") as out:
        for c, s, e, n, sc in sorted(recs):
            out.write(f"{c}\t{s}\t{e}\t{n}\t{sc}\t.\n")
print(f"{label}: {len(gained)} gained, {len(lost)} lost", file=sys.stderr)
