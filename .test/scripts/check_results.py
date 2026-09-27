#!/usr/bin/env python3
"""Regression checks on the .test outputs of `pixi run test-all`.

The simulated design (simulate_cutrun.py) has known answers:
  * TF peaks sit on the simulated TF sites, and TF libraries are mostly
    sub-nucleosomal;
  * IgG hotspots are present in every library; the IgG gate removes the peaks
    that fall on them (the IgG is 3x enriched there);
  * tfA_stim has ~1.4x more target signal per cell than tfA_ctrl, so its
    greenlist and spike-in size factors (per mapped read) are ~0.7x those of
    tfA_ctrl, and the two agree;
  * under depth normalization the unchanged sites look "lost" in stim vs
    ctrl; greenlist normalization removes most of that artefact, and
    normcheck flags the contrast.
Exit status 1 on any failure.
"""

import csv
import os
import statistics
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(HERE, "results")
REF = os.path.join(HERE, "data", "ref")
failures = []


def check(ok, msg):
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def bed(path):
    with open(path) as fh:
        return [(f[0], int(f[1]), int(f[2])) for f in (l.split("\t") for l in fh) if len(f) >= 3]


def overlap_fraction(a, b):
    if not a:
        return 0.0
    hit = sum(any(c == c2 and s < e2 and s2 < e for c2, s2, e2 in b) for c, s, e in a)
    return hit / len(a)


# 1. QC summary --------------------------------------------------------------
qc = {r["Sample"]: r for r in rows(os.path.join(R, "qc/qc_summary.tsv"))}
tf = [float(qc[s]["TF_Fraction"]) for s in qc if s.startswith("tfA")]
check(min(tf) > 0.5, f"TF libraries are mostly sub-nucleosomal (min TF_Fraction {min(tf):.2f})")
check(qc["k27_ctrl_R1"]["TF_Fraction"] == "NA", "TF fraction not scored for the histone library")
check(qc["IgG_R1"]["Role"] == "control", "IgG libraries are controls")

# 2. Peaks on truth; IgG control used ------------------------------------------
truth_tf = bed(os.path.join(REF, "truth_tf_sites.bed"))
cons = bed(os.path.join(R, "peaks/consensus/tfA.consensus.bed"))
prec = overlap_fraction(cons, truth_tf)
check(len(cons) > 100 and prec > 0.8,
      f"tfA consensus: {len(cons)} peaks, {100 * prec:.0f}% on simulated TF sites")
def macs2_log(name):
    with open(os.path.join(HERE, f"logs/peaks/macs2/merged/{name}.log")) as fh:
        return fh.read()


check("IgG.bam" in macs2_log("tfA_ctrl"), "MACS2 tfA_ctrl called against its IgG control")
check("# control file = None" in macs2_log("k27_ctrl"), "MACS2 k27_ctrl (no control) called without one")

# 3. IgG gate removes hotspot peaks ------------------------------------------
# k27_ctrl has no control, so its raw peaks include the IgG hotspots; the
# gate (against the pooled IgG) must remove them and keep the domains.
hot = bed(os.path.join(REF, "truth_hotspots.bed"))
dom = bed(os.path.join(REF, "truth_domains.bed"))
raw = bed(os.path.join(R, "peaks/macs2/merged/k27_ctrl.raw.bed"))
final = bed(os.path.join(R, "peaks/macs2/merged/k27_ctrl.bed"))
raw_hot = sum(overlap_fraction([p], hot) > 0 and overlap_fraction([p], dom) == 0 for p in raw)
fin_hot = sum(overlap_fraction([p], hot) > 0 and overlap_fraction([p], dom) == 0 for p in final)
check(raw_hot >= 5 and fin_hot == 0,
      f"IgG gate: hotspot-only peaks {raw_hot} raw -> {fin_hot} kept ({len(raw)} -> {len(final)} peaks)")
check(overlap_fraction(final, dom) > 0.9, "IgG gate keeps the histone domains")
igg = {r["Sample"]: r for r in rows(os.path.join(R, "qc/igg_enrichment.tsv"))}
check(float(igg["tfA_ctrl"]["Median_Fold_vs_IgG"]) > 2, "tfA peaks are enriched over IgG")

# 4. Greenlist and spike-in size factors agree on the global shift ------------
def per_read_ratio(method):
    sf = {r["SampleID"]: r for r in rows(os.path.join(R, f"normalization/{method}/size_factors.tsv"))}
    norm = {s: float(r["size_factor"]) / float(r["reads"]) for s, r in sf.items()}
    stim = statistics.mean(v for s, v in norm.items() if s.startswith("tfA_stim"))
    ctrl = statistics.mean(v for s, v in norm.items() if s.startswith("tfA_ctrl"))
    return stim / ctrl


g, s = per_read_ratio("greenlist"), per_read_ratio("spikein")
check(g < 0.85, f"greenlist size factors see the global increase (stim/ctrl per read {g:.2f})")
check(s < 0.85, f"spike-in size factors see the global increase (stim/ctrl per read {s:.2f})")
check(abs(g - s) / s < 0.25, f"greenlist and spike-in agree ({g:.2f} vs {s:.2f})")

# 5. Differential binding: depth artefact removed by greenlist ----------------
depth = rows(os.path.join(R, "diff/tfA/depth/summary.tsv"))[0]
green = rows(os.path.join(R, "diff/tfA/greenlist/summary.tsv"))[0]
check(int(green["Lost"]) < int(depth["Lost"]),
      f"lost peaks: depth {depth['Lost']} -> greenlist {green['Lost']}")
check(int(green["Gained"]) >= int(depth["Gained"]),
      f"gained peaks: depth {depth['Gained']} -> greenlist {green['Gained']}")
verdict = rows(os.path.join(R, "diff/tfA/normcheck/norm_verdict.tsv"))[0]
check(verdict["normalization_sensitive"] == "TRUE", "normcheck flags stim_vs_ctrl as normalization-sensitive")

# 6. Downstream outputs exist ------------------------------------------------
for p in ("heatmaps/tfA/peaks_heatmap.png", "heatmaps/contrasts/stim_vs_ctrl_heatmap.png",
          "annotate/annotation_summary.tsv", "motifs/tfA/knownResults.txt",
          "qc/multiqc/multiqc_report.html"):
    check(os.path.getsize(os.path.join(R, p)) > 0, f"{p} written")

print(f"\n{len(failures)} failure(s)")
sys.exit(1 if failures else 0)
