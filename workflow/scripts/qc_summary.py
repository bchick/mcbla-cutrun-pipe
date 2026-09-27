"""Per-library QC summary with PASS / WARN / FAIL flags.

Metrics (NA when not available; NA is never flagged):
  fragments          read pairs in the final BAM (flagstat)
  aligned_pct        bowtie2 overall alignment rate (fastq mode)
  mito_pct           mito reads / aligned reads (fastq mode)
  dup_pct            duplicate fraction (samtools markdup, fastq mode)
  frip               fragments in the group's merged peaks (targets)
  tf_fraction        fragments <= fragments.tf_max_size (TF targets only)
  igg_pass_fraction  fraction of the library's peaks >= igg.min_fold over IgG
  spikein_pct        spike-in / host fragments (reported, not flagged)

qc.thresholds gives [pass, warn]. Higher is better except mito_pct and
dup_pct, where <= pass is PASS and <= warn is WARN. Overall is the worst flag.
Flags are reported, never enforced.
"""

import csv
import re
import sys

sm = snakemake  # noqa: F821  (injected by Snakemake)
sys.stderr = open(sm.log[0], "w")

METRICS = ["fragments", "aligned_pct", "mito_pct", "dup_pct", "frip",
           "tf_fraction", "igg_pass_fraction", "spikein_pct"]
LOWER_IS_BETTER = {"mito_pct", "dup_pct"}
COLUMNS = {
    "fragments": "Fragments",
    "aligned_pct": "Aligned_Pct",
    "mito_pct": "Mito_Pct",
    "dup_pct": "Dup_Pct",
    "frip": "FRiP",
    "tf_fraction": "TF_Fraction",
    "igg_pass_fraction": "IgG_Pass_Fraction",
    "spikein_pct": "Spikein_Pct",
}
RANK = {"PASS": 0, "WARN": 1, "FAIL": 2}
thresholds = sm.params.thresholds
libs = list(sm.params.libs)
role = dict(zip(libs, sm.params.roles))
ttype = dict(zip(libs, sm.params.types))
tf_max = int(sm.params.tf_max)


def read_rows(path):
    with open(path) as fh:
        return {row["Sample"]: row for row in csv.DictReader(fh, delimiter="\t")}


def number(x):
    try:
        return float(str(x).rstrip("%"))
    except (TypeError, ValueError):
        return None


values = {lib: dict.fromkeys(METRICS) for lib in libs}

for lib, path in zip(libs, sm.input.flagstat):
    with open(path) as fh:
        m = re.search(r"^(\d+) \+ \d+ paired in sequencing", fh.read(), flags=re.MULTILINE)
    values[lib]["fragments"] = int(m.group(1)) // 2 if m else None

for lib, path in zip(libs, sm.input.frag):
    n = small = 0
    with open(path) as fh:
        for line in fh:
            size, count = (int(x) for x in line.split())
            n += count
            small += count if size <= tf_max else 0
    if n and role[lib] == "target" and ttype[lib] == "tf":
        values[lib]["tf_fraction"] = small / n

if sm.input.get("report"):
    report = read_rows(sm.input.report)
    for lib in libs:
        row = report.get(lib, {})
        values[lib]["aligned_pct"] = number(row.get("Aligned_Pct"))
        values[lib]["mito_pct"] = number(row.get("Mito_Pct"))
        values[lib]["dup_pct"] = number(row.get("Dup_Pct"))

for path in sm.input.frip:
    for lib, row in read_rows(path).items():
        values[lib]["frip"] = number(row["FRiP"])

for path in sm.input.get("igg", []):
    for lib, row in read_rows(path).items():
        values[lib]["igg_pass_fraction"] = number(row["IgG_Pass_Fraction"])

if sm.input.get("spikein"):
    for lib, row in read_rows(sm.input.spikein).items():
        values[lib]["spikein_pct"] = number(row["Spikein_Pct"])


def flag(metric, value):
    if value is None or metric not in thresholds:
        return "NA"
    pass_at, warn_at = thresholds[metric]
    if metric in LOWER_IS_BETTER:
        if value <= pass_at:
            return "PASS"
        return "WARN" if value <= warn_at else "FAIL"
    if value >= pass_at:
        return "PASS"
    return "WARN" if value >= warn_at else "FAIL"


def fmt(metric, value):
    if value is None:
        return "NA"
    if metric == "fragments":
        return str(int(value))
    if metric.endswith("_pct"):
        return f"{value:.2f}"
    return f"{value:.4g}"


header = ["Sample", "Group", "Role"]
for m in METRICS:
    header += [COLUMNS[m], f"{COLUMNS[m]}_flag"]
header += ["Overall", "Not_Passing"]

with open(sm.output[0], "w") as out:
    out.write("\t".join(header) + "\n")
    for lib, grp in zip(libs, sm.params.groups):
        row = [lib, grp, role[lib]]
        flags = {}
        for m in METRICS:
            flags[m] = flag(m, values[lib][m])
            row += [fmt(m, values[lib][m]), flags[m]]
        scored = [f for f in flags.values() if f != "NA"]
        overall = max(scored, key=RANK.get) if scored else "NA"
        not_passing = [f"{COLUMNS[m]}:{f}" for m, f in flags.items() if f in ("WARN", "FAIL")]
        row += [overall, ";".join(not_passing) or "-"]
        out.write("\t".join(row) + "\n")
        if overall == "FAIL":
            print(f"QC FAIL {lib}: {';'.join(not_passing)}", file=sys.stderr)
