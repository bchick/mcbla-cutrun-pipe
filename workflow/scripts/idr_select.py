"""Select the reproducible (IDR) peak set for one group.

With two replicates there is exactly one pair and this is a copy. With more
than two replicates every true-replicate pair was run through IDR and the pair
with the most peaks passing the threshold is kept (ENCODE ATAC "conservative"
set, without pseudo-replicates). The per-group row feeds idr_summary.tsv
with the same columns as the lab's 2.2_atac_idr.sh plus N_Reps, Selected_Pair,
All_Pairs and Status.
"""

import os
import shutil
import sys

sm = snakemake  # noqa: F821
sys.stderr = open(sm.log[0], "w")


def nlines(path):
    if not os.path.exists(path):
        return 0
    with open(path) as fh:
        return sum(1 for _ in fh)


pairs = []
for f in sm.input.pairs:
    status_file = f.replace(".idr.narrowPeak", ".status")
    status = open(status_file).read().strip() if os.path.exists(status_file) else "OK"
    name = os.path.basename(f).replace(".idr.narrowPeak", "")
    pairs.append((name, f, nlines(f), status))

ok = [p for p in pairs if p[3] == "OK"]
rep_counts = {os.path.basename(r).replace(".raw.bed", ""): nlines(r)
              for r in sm.input.reps}
oracle = nlines(sm.input.oracle)

if ok:
    best = max(ok, key=lambda p: p[2])
    shutil.copyfile(best[1], sm.output.peaks)
    a, b = best[0].split("__")
    r1, r2 = rep_counts.get(a, 0), rep_counts.get(b, 0)
    n_idr = best[2]
    avg = (r1 + r2) / 2
    reprod = f"{100 * n_idr / avg:.1f}%" if avg > 0 else "NA"
    status = "OK" if len(ok) == len(pairs) else "PARTIAL"
    selected = best[0]
else:
    open(sm.output.peaks, "w").close()
    a = b = selected = "NA"
    r1 = r2 = n_idr = 0
    reprod = "NA"
    status = "FAILED"
    print(f"WARNING: all IDR pairs failed for {sm.params.condition}", file=sys.stderr)

if len(sm.params.libs) > 2:
    if sm.params.strategy == "first_two":
        status += ";first_two"
    else:
        status += ";encode_max_pair"

all_pairs = ";".join(f"{p[0]}={p[2] if p[3] == 'OK' else 'FAILED'}" for p in pairs)
with open(sm.output.stats, "w") as out:
    out.write("\t".join(str(x) for x in [
        sm.params.condition, len(sm.params.libs), r1, r2, oracle, n_idr,
        reprod, selected, all_pairs, status]) + "\n")
