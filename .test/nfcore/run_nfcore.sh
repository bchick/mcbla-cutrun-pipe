#!/usr/bin/env bash
# run_nfcore.sh -- run nf-core/cutandrun on the synthetic .test dataset.
#
# Writes an nf-core samplesheet from .test/config/samples.tsv (absolute FASTQ
# paths; rows sharing group + replicate are runs, merged as in our pipeline)
# and runs the pinned nf-core revision with params.yaml + nfcore.config.
# Everything lands in .test/nfcore/run/ (gitignored); reruns use -resume.
# compare_nfcore.py then scores both pipelines against the simulated truth.
#
# nf-core/cutandrun 3.2.2 cannot mix controlled and uncontrolled targets: with
# use_control on, a row with an empty control is paired with no IgG and gets
# no peaks. Rows without a control are therefore given the first IgG group of
# the samplesheet (k27_ctrl in .test); the report notes this.
#
# Usage (repo root):  pixi run test-nfcore-run
#   NFCORE_PROFILE=singularity pixi run test-nfcore-run   (default: docker)
#   NFCORE_REVISION=3.2.2                                 (pinned default)
#
# Nextflow is pinned to 24.04.4 (NXF_VER; downloaded on first use): 3.2.2's
# Trim Galore module uses `emit: html optional true`, which Nextflow 25.x
# rejects ("Cannot invoke method optional() on null object").
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd)
HERE=$REPO/.test/nfcore
RUN=$HERE/run
PROFILE=${NFCORE_PROFILE:-docker}
REVISION=${NFCORE_REVISION:-3.2.2}
DATA=$REPO/.test/data
export NXF_VER=${NXF_VER:-24.04.4}

command -v nextflow >/dev/null || { echo "nextflow not on PATH (see .test/README.md)" >&2; exit 1; }
[[ -s $DATA/ref/genome.fa ]] || { echo "no test data; run: pixi run build-test" >&2; exit 1; }

mkdir -p "$RUN"
awk -F'\t' -v d="$DATA/fastq" '
    NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    {n++; g[n] = $c["group"]; r[n] = $c["replicate"]; f1[n] = $c["fastq_1"]; f2[n] = $c["fastq_2"]
     ctl[n] = $c["control"]; if (ctl[n] != "") used[ctl[n]] = 1}
    END {
        for (i = 1; i <= n; i++) if (g[i] in used) {igg = g[i]; break}
        OFS = ","; print "group", "replicate", "fastq_1", "fastq_2", "control"
        for (i = 1; i <= n; i++) {
            x = ctl[i]; if (x == "" && !(g[i] in used)) x = igg
            print g[i], r[i], d "/" f1[i], d "/" f2[i], x
        }
    }' "$REPO/.test/config/samples.tsv" > "$RUN/samplesheet.csv"

# 3.2.2's blacklist step (BEDTOOLS_INTERSECT, meta.id "blacklist") refuses an
# input named blacklist.bed ("Input and output names are the same").
ln -sf "$DATA/ref/blacklist.bed" "$RUN/test_blacklist.bed"

cd "$RUN"
nextflow run nf-core/cutandrun -r "$REVISION" -profile "$PROFILE" \
    -params-file "$HERE/params.yaml" -c "$HERE/nfcore.config" \
    --input "$RUN/samplesheet.csv" --outdir "$RUN/results" \
    --fasta "$DATA/ref/genome.fa" --gtf "$DATA/ref/genes.gtf" \
    --blacklist "$RUN/test_blacklist.bed" --spikein_fasta "$DATA/ref/spikein.fa" \
    -resume
