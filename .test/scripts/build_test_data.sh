#!/usr/bin/env bash
# Build the synthetic .test dataset (idempotent; outputs are gitignored).
#
#   bash .test/scripts/build_test_data.sh [--fasta hg38.fa] [--pairs N] [--force]
#
# With --fasta, chr22 and chrM are read from a local hg38/GRCh38 FASTA
# (plain FASTA with .fai, or .gz). Without it, chr22 and chrM are downloaded
# from UCSC (~12 MB). Only Python 3 (standard library) and curl are needed;
# the bowtie2 index is built by the pipeline itself (reference.bowtie2_index: "").
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$HERE/data"
FASTA=""
PAIRS=60000
FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta) FASTA="$2"; shift 2 ;;
        --pairs) PAIRS="$2"; shift 2 ;;
        --force) FORCE=1; shift ;;
        -h|--help) sed -n '2,10p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -s "$OUT/ref/genome.fa" && -s "$OUT/fastq/IgG_r2_R2.fastq.gz" && $FORCE -eq 0 ]]; then
    echo "Test data already present in $OUT (use --force to rebuild)"
    exit 0
fi
rm -rf "$OUT"; mkdir -p "$OUT/download"

CHRM_ARG=()
if [[ -z "$FASTA" ]]; then
    UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes
    for c in chr22 chrM; do
        echo "Downloading $c from UCSC..."
        curl -fsSL -o "$OUT/download/$c.fa.gz" "$UCSC/$c.fa.gz"
    done
    FASTA="$OUT/download/chr22.fa.gz"
    CHRM_ARG=(--chrm-fasta "$OUT/download/chrM.fa.gz")
fi

echo "Simulating reads ($PAIRS pairs per run)..."
python3 "$HERE/scripts/simulate_cutrun.py" --fasta "$FASTA" "${CHRM_ARG[@]}" \
    --outdir "$OUT" --pairs "$PAIRS"
rm -rf "$OUT/download"
du -sh "$OUT"/*
echo "Done. Run the test from the repo root with: pixi run test"
