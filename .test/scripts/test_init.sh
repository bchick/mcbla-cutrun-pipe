#!/usr/bin/env bash
# test_init.sh -- check workflow/scripts/init_project.py against a fixture manifest.
#
# Builds a throwaway manifest whose "resources" are empty placeholder files,
# runs init non-interactively for each genome x IgG x spike-in x blacklist
# choice, and dry-runs the Snakemake DAG of every generated project. Also checks
# that init refuses to guess when answers are missing.
#
# Usage (repo root):  pixi run test-init
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd)
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT

# --- placeholder resources ---------------------------------------------------
R=$T/resource
mkdir -p "$R"
for g in hs mm; do
    touch "$R/$g.fa" "$R/$g.fa.fai" "$R/$g.gtf" "$R/$g.blacklist.bed" "$R/$g.greenlist.bed"
    for s in 1 2 3 4 rev.1 rev.2; do touch "$R/$g.$s.bt2" "$R/${g}_ecoli.$s.bt2"; done
done
cat > "$T/manifest.yaml" <<EOF
manifest_version: 1
updated: test
genomes:
  hg38:
    description: "fixture human"
    aliases: [GRCh38]
    fasta: $R/hs.fa
    bowtie2_index: $R/hs
    gtf: $R/hs.gtf
    effective_genome_size: 2913022398
    macs2_gsize: hs
    mito_chrom: chrM
    blacklists:
      encode_v2: {path: $R/hs.blacklist.bed, default: true}
    greenlist: {path: $R/hs.greenlist.bed}
    spikeins:
      ecoli: {combined_bowtie2_index: $R/hs_ecoli, contig_prefix: "Ecoli_"}
      dm6: {combined_bowtie2_index: $R/nonexistent, contig_prefix: "dm6_", status: unverified}
  mm10:
    description: "fixture mouse, no greenlist"
    fasta: $R/mm.fa
    bowtie2_index: $R/mm
    gtf: $R/mm.gtf
    effective_genome_size: 2652783500
    macs2_gsize: mm
    blacklists:
      encode_v2: {path: $R/mm.blacklist.bed, default: true}
    spikeins:
      ecoli: {combined_bowtie2_index: $R/mm_ecoli, contig_prefix: "Ecoli_"}
EOF

init() { python "$REPO/workflow/scripts/init_project.py" --manifest "$T/manifest.yaml" "$@" < /dev/null; }
fail() { echo "FAIL: $*" >&2; exit 1; }

# --- refusals ------------------------------------------------------------------
set +e
init --dir "$T/x" --genome hg38 > /dev/null 2>&1; [[ $? == 2 ]] || fail "missing answers should exit 2"
init --dir "$T/x" --genome hg38 --igg yes --spikein dm6 --blacklist none > /dev/null 2>&1
[[ $? == 1 ]] || fail "unverified spike-in should be refused"
init --dir "$T/x" --genome hg19 --igg yes --spikein none --blacklist none > /dev/null 2>&1
[[ $? == 1 ]] || fail "unknown genome should be refused"
set -e
init --list > /dev/null
init --list --json | python -c "import json,sys; d=json.load(sys.stdin); assert 'dm6' not in d['hg38']['spikeins']"

# --- every combination dry-runs -------------------------------------------------
n=0
for genome in hg38 mm10; do
    for igg in yes no; do
        for spike in ecoli none; do
            for bl in encode_v2 none; do
                d=$T/proj_${genome}_${igg}_${spike}_${bl}
                init --dir "$d" --genome "$genome" --igg "$igg" --spikein "$spike" --blacklist "$bl" > /dev/null
                mkdir -p "$d/fastq"
                if [[ $igg == no ]]; then  # drop the IgG rows and controls from the example sheet
                    awk -F'\t' 'BEGIN{OFS="\t"} NR == 1 || $6 != "IgG" {if (NR > 1) $5 = ""; print}' \
                        "$d/samples.tsv" > "$d/s" && mv "$d/s" "$d/samples.tsv"
                fi
                for f in $(tail -n +2 "$d/samples.tsv" | cut -f3,4); do touch "$d/fastq/$f"; done
                snakemake -n -s "$REPO/workflow/Snakefile" --directory "$d" \
                    --configfile "$d/project.yaml" > "$d.log" 2>&1 \
                    || { cat "$d.log"; fail "dry run $genome igg=$igg spikein=$spike blacklist=$bl"; }
                grep -q "spikein_counts" "$d.log" && [[ $spike == none ]] && fail "spike-in jobs without spike-in"
                grep -q "igg_gate" "$d.log" && [[ $igg == no ]] && fail "IgG jobs without IgG"
                n=$((n + 1))
            done
        done
    done
done
echo "test_init: refusals OK, $n generated projects dry-run OK"
