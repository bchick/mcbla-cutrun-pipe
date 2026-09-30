#!/usr/bin/env python3
"""Head-to-head of this pipeline and nf-core/cutandrun on the synthetic .test data.

Run from the repo root after `pixi run test-all` and `pixi run test-nfcore-run`:

    python .test/nfcore/compare_nfcore.py

Both pipelines are scored against the simulated truth (simulate_cutrun.py):
TF sites (truth_tf_sites.bed), H3K27-like domains (truth_domains.bed), IgG
hotspots (truth_hotspots.bed) and 8% PCR duplicates. Two kinds of check:

  gate    fails the test (exit 1) when this pipeline is meaningfully worse
          than nf-core against the truth:
            * tfA consensus peaks: F1 or precision more than 0.05 below
              nf-core's;
            * k27 domain recall more than 0.05 below nf-core's;
            * k27 peaks on IgG hotspots only (artefacts): more than nf-core;
            * per-library duplicate rate: further from the simulated value
              than nf-core's, by more than 0.02.
  report  written to the report only (peak counts, Jaccard, FRiP, spike-in).

nf-core/cutandrun 3.2.2 cannot leave one target uncontrolled while others use
IgG, so run_nfcore.sh gives k27_ctrl the IgG control there. Our pipeline calls
k27_ctrl without a control and removes hotspots with the IgG gate instead.

Writes .test/results/nfcore_compare/{metrics.tsv,report.md}.
"""

import bisect
import csv
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OURS = os.path.join(HERE, "results")
NF = os.path.join(HERE, "nfcore", "run", "results")
REF = os.path.join(HERE, "data", "ref")
OUT = os.path.join(OURS, "nfcore_compare")

SIM_DUP = 0.08
TOL_PEAKS = 0.05
TOL_QC = 0.02

metrics = []  # (section, metric, ours, nfcore, kind, status)
failures = []


def fmt(x):
    return "NA" if x is None else (f"{x:.3f}" if isinstance(x, float) else str(x))


def record(section, metric, ours, nfcore, gate=None):
    """gate: None (report) or a bool (True = pass)."""
    kind = "report" if gate is None else "gate"
    status = "-" if gate is None else ("PASS" if gate else "FAIL")
    metrics.append((section, metric, ours, nfcore, kind, status))
    print(f"{status:4} {section}: {metric}  ours={fmt(ours)}  nf-core={fmt(nfcore)}")
    if gate is False:
        failures.append(f"{section}: {metric}")


def need(pattern, what):
    hits = sorted(glob.glob(pattern))
    if not hits:
        sys.exit(f"missing {what}: {pattern}\n(run `pixi run test-all` and `pixi run test-nfcore-run` first)")
    return hits


# --- intervals ----------------------------------------------------------------
def bed(path):
    with open(path) as fh:
        return [(f[0], int(f[1]), int(f[2])) for f in (l.split("\t") for l in fh)
                if len(f) >= 3 and not f[0].startswith(("#", "track"))]


def merged(ivs):
    """Disjoint sorted intervals per chromosome."""
    by = {}
    for c, s, e in sorted(ivs):
        m = by.setdefault(c, [])
        if m and s <= m[-1][1]:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    return {c: ([s for s, _ in m], [e for _, e in m]) for c, m in by.items()}


def hits(a, b):
    """For each interval in a, whether it overlaps any interval in b."""
    mb = merged(b)
    out = []
    for c, s, e in a:
        ok = False
        if c in mb:
            starts, ends = mb[c]
            i = bisect.bisect_left(starts, e) - 1
            ok = i >= 0 and ends[i] > s
        out.append(ok)
    return out


def hit_fraction(a, b):
    return sum(hits(a, b)) / len(a) if a else 0.0


def jaccard(a, b):
    """Base-pair Jaccard of two interval sets."""
    def total(m):
        return sum(e - s for st, en in m.values() for s, e in zip(st, en))
    ma, mb = merged(a), merged(b)
    inter = 0
    for c in ma.keys() & mb.keys():
        (sa, ea), (sb, eb) = ma[c], mb[c]
        i = j = 0
        while i < len(sa) and j < len(sb):
            inter += max(0, min(ea[i], eb[j]) - max(sa[i], sb[j]))
            if ea[i] < eb[j]:
                i += 1
            else:
                j += 1
    union = total(ma) + total(mb) - inter
    return inter / union if union else 0.0


def prf(peaks, truth):
    p, r = hit_fraction(peaks, truth), hit_fraction(truth, peaks)
    return p, r, (2 * p * r / (p + r) if p + r else 0.0)


def hotspot_only(peaks, hot, dom):
    """Peaks on an IgG hotspot and on no true domain."""
    return sum(h and not d for h, d in zip(hits(peaks, hot), hits(peaks, dom)))


def nf_consensus(group):
    """nf-core per-group consensus (replicate-threshold filtered when written)."""
    pat = os.path.join(NF, f"03_peak_calling/05_consensus_peaks/{group}.macs2.consensus.peak_counts*.bed")
    found = need(pat, f"nf-core consensus peaks for {group}")
    filtered = [f for f in found if f.endswith(".awk.bed")]
    return bed((filtered or found)[0])


# --- peaks --------------------------------------------------------------------
tf = bed(os.path.join(REF, "truth_tf_sites.bed"))
dom = bed(os.path.join(REF, "truth_domains.bed"))
hot = bed(os.path.join(REF, "truth_hotspots.bed"))

# tfA: ours is one consensus per target; nf-core's is per group, so its
# tfA_ctrl and tfA_stim consensus sets are pooled.
ours_tf = bed(need(os.path.join(OURS, "peaks/consensus/tfA.consensus.bed"), "our tfA consensus")[0])
nf_tf = nf_consensus("tfA_ctrl") + nf_consensus("tfA_stim")
po, ro, fo = prf(ours_tf, tf)
pn, rn, fn = prf(nf_tf, tf)
record("tfA consensus", "peaks", len(ours_tf), sum(len(st) for st, _ in merged(nf_tf).values()))
record("tfA consensus", "precision vs TF sites", po, pn, po >= pn - TOL_PEAKS)
record("tfA consensus", "recall vs TF sites", ro, rn)
record("tfA consensus", "F1 vs TF sites", fo, fn, fo >= fn - TOL_PEAKS)
record("tfA consensus", "bp Jaccard ours vs nf-core", jaccard(ours_tf, nf_tf), None)

# k27: broad domains; nf-core calls it against IgG, ours without + IgG gate.
ours_k27 = bed(need(os.path.join(OURS, "peaks/consensus/k27.consensus.bed"), "our k27 consensus")[0])
ours_k27_final = bed(need(os.path.join(OURS, "peaks/macs2/merged/k27_ctrl.bed"), "our gated k27 peaks")[0])
nf_k27 = nf_consensus("k27_ctrl")
ro, rn = hit_fraction(dom, ours_k27), hit_fraction(dom, nf_k27)
record("k27 consensus", "peaks", len(ours_k27), len(nf_k27))
record("k27 consensus", "domain recall", ro, rn, ro >= rn - TOL_PEAKS)
record("k27 consensus", "precision vs domains", hit_fraction(ours_k27, dom), hit_fraction(nf_k27, dom))
ho, hn = hotspot_only(ours_k27_final, hot, dom), hotspot_only(nf_k27, hot, dom)
record("k27 peaks", "hotspot-only peaks kept (IgG artefacts)", ho, hn, ho <= hn)
record("k27 consensus", "bp Jaccard ours vs nf-core", jaccard(ours_k27, nf_k27), None)


# --- per-library QC -----------------------------------------------------------
def ours_qc():
    with open(os.path.join(OURS, "qc/qc_summary.tsv")) as fh:
        return {r["Sample"]: r for r in csv.DictReader(fh, delimiter="\t")}


def nf_dup(lib):
    """Duplicate fraction from nf-core's post-MarkDuplicates flagstat."""
    path = need(os.path.join(NF, f"02_alignment/*/target/markdup/{lib}.flagstat"),
                f"nf-core markdup flagstat for {lib}")[0]
    with open(path) as fh:
        t = fh.read()
    n = lambda k: int(re.search(rf"^(\d+) \+ \d+ {k}", t, re.M).group(1))
    return n("duplicates") / n("in total")


def concordant_pairs(path):
    with open(path) as fh:
        t = fh.read()
    return sum(int(x) for x in re.findall(r"^\s*(\d+) \([\d.]+%\) aligned concordantly (?:exactly 1|>1) time", t, re.M))


def nf_spike_ratio(lib):
    """Spike-in / host concordant pairs from nf-core's bowtie2 logs."""
    sp = glob.glob(os.path.join(NF, f"02_alignment/*/spikein/log/{lib}.spikein.bowtie2.log"))
    tg = glob.glob(os.path.join(NF, f"02_alignment/*/target/log/{lib}.bowtie2.log"))
    if not sp or not tg:
        return None
    return concordant_pairs(sp[0]) / concordant_pairs(tg[0])


def frip(lib):
    p = os.path.join(OURS, f"qc/frip/{lib}.frip.tsv")
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return float(next(csv.DictReader(fh, delimiter="\t"))["FRiP"])


qc = ours_qc()
for lib, row in sorted(qc.items()):
    do, dn = float(row["Dup_Pct"]) / 100, nf_dup(lib)
    record(lib, f"duplicate rate (simulated {SIM_DUP})", do, dn,
           abs(do - SIM_DUP) <= abs(dn - SIM_DUP) + TOL_QC)
    if row["Role"] != "control":
        record(lib, "FRiP", frip(lib), None)

# Spike-in: tfA_stim has ~1.4x the target signal per cell of tfA_ctrl, so its
# spike-in share per host read drops (stim/ctrl ~0.7). Report only.
with open(os.path.join(OURS, "qc/spikein_summary.tsv")) as fh:
    sp = {r["Sample"]: int(r["Spikein_Fragments"]) / int(r["Host_Fragments"])
          for r in csv.DictReader(fh, delimiter="\t")}


def ratio(get):
    vals = {s: get(s) for s in qc if s.startswith("tfA")}
    if any(v is None for v in vals.values()):
        return None
    mean = lambda pre: sum(v for s, v in vals.items() if s.startswith(pre)) / sum(s.startswith(pre) for s in vals)
    return mean("tfA_stim") / mean("tfA_ctrl")


record("spike-in", "tfA stim/ctrl spike-in per host read (expected < 0.85)", ratio(sp.get), ratio(nf_spike_ratio))


# --- write --------------------------------------------------------------------
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "metrics.tsv"), "w") as fh:
    fh.write("section\tmetric\tours\tnfcore\tkind\tstatus\n")
    for m in metrics:
        fh.write("\t".join(fmt(x) for x in m) + "\n")
with open(os.path.join(OUT, "report.md"), "w") as fh:
    fh.write("# mcbla-cutrun-pipe vs nf-core/cutandrun on the .test dataset\n\n")
    fh.write("Truth: `.test/data/ref/truth_{tf_sites,domains,hotspots}.bed`, 8% PCR duplicates.\n")
    fh.write("Gates fail the test; report rows are for reading only. nf-core calls k27_ctrl against "
             "IgG (3.2.2 cannot leave one target uncontrolled); ours calls it without a control and "
             "applies the IgG gate.\n\n")
    fh.write(f"**{len(failures)} gate failure(s).**\n\n")
    fh.write("| section | metric | ours | nf-core | kind | status |\n|---|---|---|---|---|---|\n")
    for m in metrics:
        fh.write("| " + " | ".join(fmt(x) for x in m) + " |\n")
print(f"\nwrote {OUT}/metrics.tsv and report.md")
print(f"{len(failures)} gate failure(s)")
sys.exit(1 if failures else 0)
