"""Greenlist or spike-in size factors per normalization group.

greenlist: DESeq2 median-of-ratios on the greenlist count matrix, computed
           separately within each normalization group (target, or target x
           batch). Same formula as DESeq2::estimateSizeFactorsForMatrix: rows
           with a zero in any library of the group are dropped, the factor is
           the median ratio of a library to the row geometric mean.
spikein:   spike-in fragments / geometric mean of spike-in fragments in the group.

Size factors are computed inside each group because factors estimated over
a larger set than the one being compared differ (mcf7 SD51: r = 0.85 between
project-wide and subset factors). A group with one library gets 1.

Columns: SampleID, group, target, norm_group, size_factor, reads, scale_factor
where reads are mapped reads (flagstat) and
scale_factor = 1e6 / (size_factor * geometric mean of reads in the group),
the bamCoverage --scaleFactor giving CPM-equivalent tracks.
"""

import csv
import math
import re
import statistics
import sys
from collections import defaultdict

sm = snakemake  # noqa: F821
sys.stderr = open(sm.log[0], "w")

libs = list(sm.params.libs)
meta = dict(zip(libs, zip(sm.params.groups, sm.params.targets, sm.params.norm_groups)))


def geomean(xs):
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


reads = {}
for lib, path in zip(libs, sm.input.flagstat):
    with open(path) as fh:
        m = re.search(r"^(\d+) \+ \d+ mapped", fh.read(), flags=re.MULTILINE)
    reads[lib] = int(m.group(1)) if m else 0

members = defaultdict(list)
for lib in libs:
    members[meta[lib][2]].append(lib)

sf = {}
if sm.params.method == "greenlist":
    with open(sm.input.counts) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        labels = [h.strip("#'") for h in header[3:]]
        matrix = {lab: [] for lab in labels}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            for lab, v in zip(labels, f[3:]):
                matrix[lab].append(float(v) if v != "nan" else 0.0)
    missing = set(libs) - set(matrix)
    if missing:
        sys.exit(f"greenlist counts lack libraries: {sorted(missing)}")
    for grp, ls in members.items():
        if len(ls) == 1:
            print(f"WARNING: normalization group {grp} has one library; "
                  "size factor 1 (CPM)", file=sys.stderr)
            sf[ls[0]] = 1.0
            continue
        rows = list(zip(*(matrix[l] for l in ls)))
        usable = [r for r in rows if all(x > 0 for x in r)]
        print(f"{grp}: {len(usable)} of {len(rows)} greenlist regions with "
              "counts in every library", file=sys.stderr)
        if len(usable) < 50:
            print(f"WARNING: {grp}: only {len(usable)} usable greenlist regions; "
                  "size factors are noisy (is the greenlist for this genome?)",
                  file=sys.stderr)
        if not usable:
            sys.exit(f"{grp}: no greenlist region has counts in every library")
        gms = [geomean(r) for r in usable]
        for j, lib in enumerate(ls):
            sf[lib] = statistics.median(r[j] / g for r, g in zip(usable, gms))
else:
    spike = {}
    with open(sm.input.counts) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            spike[row["Sample"]] = float(row["Spikein_Fragments"])
    for grp, ls in members.items():
        zero = [l for l in ls if spike.get(l, 0) <= 0]
        if zero:
            sys.exit(f"{grp}: no spike-in fragments for {', '.join(zero)}")
        g = geomean([spike[l] for l in ls])
        for lib in ls:
            sf[lib] = spike[lib] / g

with open(sm.output[0], "w") as out:
    out.write("SampleID\tgroup\ttarget\tnorm_group\tsize_factor\treads\tscale_factor\n")
    for lib in libs:
        grp = meta[lib][2]
        gm_reads = geomean([max(reads[l], 1) for l in members[grp]])
        scale = 1e6 / (sf[lib] * gm_reads)
        out.write(f"{lib}\t{meta[lib][0]}\t{meta[lib][1]}\t{grp}\t{sf[lib]:.6f}\t"
                  f"{reads[lib]}\t{scale:.6g}\n")
        print(f"{lib}\t{grp}\tsf={sf[lib]:.4f}\treads={reads[lib]}\tscale={scale:.4g}",
              file=sys.stderr)
