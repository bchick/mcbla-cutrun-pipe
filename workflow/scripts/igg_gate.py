"""IgG enrichment gate for one peak set.

Per peak, target and IgG fragment centres are counted (multiBamSummary
BED-file --centerReads) together with a shared set of background windows.
Each library is scaled by its own background density (fragments per bp over
the windows), so the fold is free of library depth and of the target's FRiP:

    fold = ((t + 1) / d_target) / ((g + 1) / d_igg)

with t, g the peak counts and d the background densities; the peak width
cancels. This is the design of the mcf7 analysis 39 "IgG reality gate"
(enrichment over per-library background, one counting instrument for target
and IgG, threshold 2x).

Outputs
  <name>.igg.tsv          per peak: coordinates, counts, fold, hotspot overlap, pass
  <name>.igg_summary.tsv  one row: peaks, median fold, pass fraction, hotspot
                          overlap fraction, peaks kept
  <name>.bed              the peaks kept: all of them unless igg.hotspot_filter
                          is on, then fold >= igg.min_fold (and, with
                          igg.remove_hotspot_overlaps, no IgG hotspot overlap)
"""

import statistics
import sys

sm = snakemake  # noqa: F821
sys.stderr = open(sm.log[0], "w")

min_fold = float(sm.params.min_fold)


def key(chrom, start, end):
    return (chrom, int(start), int(end))


counts = {}
with open(sm.input.counts) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        counts[key(*f[:3])] = (float(f[3]), float(f[4]))

bg_t = bg_g = bg_bp = 0.0
with open(sm.input.bg) as fh:
    for line in fh:
        f = line.split("\t")
        k = key(*f[:3])
        if k in counts:
            t, g = counts[k]
            bg_t += t
            bg_g += g
            bg_bp += k[2] - k[1]

d_t = (bg_t + 1) / bg_bp if bg_bp else 1.0
d_g = (bg_g + 1) / bg_bp if bg_bp else 1.0
print(f"background: {bg_bp:.0f} bp, target {bg_t:.0f}, IgG {bg_g:.0f} fragments",
      file=sys.stderr)

hot = {}
with open(sm.input.hot) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4:
            hot[key(*f[:3])] = int(f[3]) > 0

rows, folds, kept = [], [], []
n_hot = 0
with open(sm.input.peaks) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        k = key(*f[:3])
        t, g = counts.get(k, (0.0, 0.0))
        fold = ((t + 1) / d_t) / ((g + 1) / d_g)
        on_hot = hot.get(k, False)
        n_hot += on_hot
        passed = fold >= min_fold
        keep = True
        if sm.params.filter:
            keep = passed and not (sm.params.remove_hotspots and on_hot)
        folds.append(fold)
        rows.append([*f[:3], f"{t:.0f}", f"{g:.0f}", f"{fold:.3f}", str(int(on_hot)),
                     str(int(passed)), str(int(keep))])
        if keep:
            kept.append(line)

with open(sm.output.table, "w") as out:
    out.write("chrom\tstart\tend\ttarget_frags\tigg_frags\tfold\tigg_hotspot\t"
              "pass\tkept\n")
    for r in rows:
        out.write("\t".join(r) + "\n")

with open(sm.output.peaks, "w") as out:
    out.writelines(kept)

n = len(rows)
med = f"{statistics.median(folds):.3f}" if folds else "NA"
pass_frac = f"{sum(f >= min_fold for f in folds) / n:.4f}" if n else "NA"
hot_frac = f"{n_hot / n:.4f}" if n else "NA"
with open(sm.output.summary, "w") as out:
    out.write("Sample\tPeaks\tMedian_Fold_vs_IgG\tIgG_Pass_Fraction\t"
              "IgG_Hotspot_Fraction\tPeaks_Kept\n")
    out.write(f"{sm.params.name}\t{n}\t{med}\t{pass_frac}\t{hot_frac}\t{len(kept)}\n")
print(f"{sm.params.name}: {n} peaks, median fold {med}, pass {pass_frac}, "
      f"kept {len(kept)}", file=sys.stderr)
