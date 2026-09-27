"""Write the DiffBind sample sheet for one target (one row per library).

Condition is the samplesheet group. Factor is the batch when config `batch`
is true (design ~Factor + Condition), otherwise the target.
"""

import csv
import sys

sm = snakemake  # noqa: F821
sys.stderr = open(sm.log[0], "w")

with open(sm.output[0], "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["SampleID", "Tissue", "Factor", "Condition", "Treatment",
                "Replicate", "bamReads", "Peaks", "PeakCaller"])
    for lib, grp, rep, factor, bam, peaks in zip(
            sm.params.libs, sm.params.groups, sm.params.replicates,
            sm.params.factors, sm.input.bams, sm.input.peaks):
        w.writerow([lib, sm.params.tissue, factor, grp, grp, rep, bam, peaks,
                    "narrow"])
