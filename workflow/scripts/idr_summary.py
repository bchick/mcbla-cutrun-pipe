"""Collect per-group IDR rows into idr_summary.tsv."""

import sys

sm = snakemake  # noqa: F821
sys.stderr = open(sm.log[0], "w")

header = ["Group", "N_Reps", "Rep1_Peaks", "Rep2_Peaks", "Oracle_Peaks",
          "IDR_Peaks", "Reprod_Rate", "Selected_Pair", "All_Pairs", "Status"]

with open(sm.output[0], "w") as out:
    out.write("\t".join(header) + "\n")
    for f in sm.input.stats:
        with open(f) as fh:
            out.write(fh.read())
    for cond, peaks in zip(sm.params.single, sm.input.get("single", [])):
        with open(peaks) as fh:
            n = sum(1 for _ in fh)
        status = ("single_rep:stringent_fallback" if sm.params.fallback == "stringent"
                  else "single_rep:skipped")
        used = n if sm.params.fallback == "stringent" else 0
        out.write("\t".join(str(x) for x in [
            cond, 1, "NA", "NA", "NA", used, "NA", "NA", "NA", status]) + "\n")
