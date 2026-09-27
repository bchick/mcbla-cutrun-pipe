#!/usr/bin/env python3
"""Simulate a tiny paired-end CUT&RUN dataset for the pipeline test.

Standard library only. Produces, under --outdir:
  ref/genome.fa        a 4 Mb window of chr22 (renamed "chr22") + chrM
  ref/spikein.fa       a random 300 kb "E. coli" contig (spike-in carry-over)
  ref/genes.gtf        synthetic transcripts whose TSSs sit on TF sites
  ref/blacklist.bed    two synthetic blacklist regions (with artefact pileups)
  ref/greenlist.bed    40 regions with the same absolute signal in every library
  ref/truth_*.bed      simulated TF sites, histone domains and IgG hotspots
  fastq/*.fastq.gz     paired-end 50 bp reads, one pair of files per run

Model. Every library is a mixture of components with an *absolute* weight
per cell: genome background, greenlist regions, IgG hotspots (sticky
regions present in IgG and in every target), blacklist pileups, chrM, the
spike-in and the target signal. Reads are drawn in proportion to the
weights, so a larger target signal dilutes everything else. The tfA_stim
group induces 30 % of TF sites 4-fold, a global increase in binding: depth
normalization hides it, greenlist and spike-in scaling recover it.

  tfA_ctrl   TF, 2 replicates      control IgG
  tfA_stim   TF, 3 replicates (replicate 3 sequenced in two runs), control IgG
  k27_ctrl   histone, broad, 1 replicate, no control (exercises the IgG gate)
  IgG        2 replicates

TF fragments are mostly sub-nucleosomal (< 120 bp); histone fragments are
mono/di-nucleosomal. Fragments shorter than the read length carry TruSeq
adapter read-through, so cutadapt has work to do. Coordinates are relative
to the extracted window, not to hg38.
"""

import argparse
import gzip
import math
import os
import random

ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
READ_LEN = 50
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# (group, replicate label, run label, kind)
RUNS = [
    ("tfA_ctrl", "r1", "r1", "tf"),
    ("tfA_ctrl", "r2", "r2", "tf"),
    ("tfA_stim", "r1", "r1", "tf_stim"),
    ("tfA_stim", "r2", "r2", "tf_stim"),
    ("tfA_stim", "r3", "r3a", "tf_stim"),
    ("tfA_stim", "r3", "r3b", "tf_stim"),
    ("k27_ctrl", "r1", "r1", "histone"),
    ("IgG", "r1", "r1", "igg"),
    ("IgG", "r2", "r2", "igg"),
]

# absolute weights per cell (arbitrary units)
W_BACKGROUND = 0.1   # ~0.1x coverage, as in real CUT&RUN
W_GREENLIST = 0.15
W_HOTSPOT = 0.10
W_BLACKLIST = 0.02
W_MITO = 0.03
W_SPIKE = 0.06
W_TF = 1.2          # all TF sites at baseline
W_HISTONE = 2.0
INDUCED_FRACTION = 0.3
INDUCTION = 4.0


def read_fasta_region(path, contig, start=None, end=None):
    """Return sequence of contig[start:end]; uses .fai when present."""
    fai = path + ".fai"
    if os.path.exists(fai) and not path.endswith(".gz"):
        with open(fai) as fh:
            for line in fh:
                name, length, offset, lb, lw = line.split("\t")[:5]
                if name == contig:
                    length, offset, lb, lw = int(length), int(offset), int(lb), int(lw)
                    break
            else:
                raise SystemExit(f"{contig} not in {fai}")
        s = 0 if start is None else start
        e = length if end is None else min(end, length)
        first = offset + (s // lb) * lw + s % lb
        last = offset + (e // lb) * lw + e % lb
        with open(path, "rb") as fh:
            fh.seek(first)
            raw = fh.read(last - first).decode()
        return raw.replace("\n", "").replace("\r", "").upper()
    opener = gzip.open if path.endswith(".gz") else open
    seq, keep = [], False
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if keep:
                    break
                keep = line[1:].split()[0] == contig
                continue
            if keep:
                seq.append(line.strip())
    if not seq:
        raise SystemExit(f"{contig} not found in {path}")
    return "".join(seq).upper()[start:end]


def revcomp(s):
    return s.translate(COMP)[::-1]


def mutate(s, rng, rate=0.001):
    out = list(s)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice("ACGT")
    return "".join(out)


def frag_len(rng, kind):
    if kind == "tf":
        mu, sd = (70, 15) if rng.random() < 0.8 else (180, 20)
    elif kind == "histone":
        mu, sd = (170, 20) if rng.random() < 0.8 else (330, 30)
    else:
        mu, sd = 200, 60
    return int(min(800, max(30, rng.gauss(mu, sd))))


def write_fasta(path, records):
    with open(path, "w") as out:
        for name, s in records:
            out.write(f">{name}\n")
            for i in range(0, len(s), 60):
                out.write(s[i:i + 60] + "\n")


class Sampler:
    """Weighted choice over (centre, weight) with a cumulative table."""

    def __init__(self, items):
        self.items = [c for c, _ in items]
        self.cum, t = [], 0.0
        for _, w in items:
            t += w
            self.cum.append(t)
        self.total = t

    def pick(self, rng):
        r = rng.random() * self.total
        lo, hi = 0, len(self.cum) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if self.cum[mid] < r:
                lo = mid + 1
            else:
                hi = mid
        return self.items[lo]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="FASTA with chr22 and chrM (plain+.fai, or .gz)")
    ap.add_argument("--chrm-fasta", default=None, help="separate FASTA holding chrM (optional)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--start", type=int, default=20_000_000)
    ap.add_argument("--length", type=int, default=4_000_000)
    ap.add_argument("--pairs", type=int, default=60_000, help="read pairs per run")
    ap.add_argument("--tf-sites", type=int, default=300)
    ap.add_argument("--domains", type=int, default=80)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()

    rng = random.Random(a.seed)
    ref = os.path.join(a.outdir, "ref")
    fq = os.path.join(a.outdir, "fastq")
    os.makedirs(ref, exist_ok=True)
    os.makedirs(fq, exist_ok=True)

    chrom = read_fasta_region(a.fasta, "chr22", a.start, a.start + a.length)
    chrm = read_fasta_region(a.chrm_fasta or a.fasta, "chrM")
    spike = "".join(rng.choice("ACGT") for _ in range(300_000))
    L = len(chrom)
    write_fasta(os.path.join(ref, "genome.fa"), [("chr22", chrom), ("chrM", chrm)])
    write_fasta(os.path.join(ref, "spikein.fa"), [("Ecoli_sim", spike)])

    def ok(c, half=1_500):
        s = chrom[c - half:c + half]
        return s.count("N") < 0.02 * len(s)

    taken = []

    def place(n, spacing, half):
        out = []
        for _ in range(200_000):
            if len(out) == n:
                break
            c = rng.randrange(20_000, L - 20_000)
            if ok(c, half) and all(abs(c - x) > spacing for x in taken):
                out.append(c)
                taken.append(c)
        else:
            raise SystemExit(f"could not place {n} sites {spacing} bp apart; "
                             "use a longer --length or fewer sites")
        return sorted(out)

    blacklist = place(2, 50_000, 2_500)
    hotspots = place(15, 5_000, 500)
    greenlist = place(40, 5_000, 500)
    tf_sites = place(a.tf_sites, 2_000, 500)
    domains = place(a.domains, 4_000, 2_000)

    tf_strength = [math.exp(rng.gauss(0, 0.6)) for _ in tf_sites]
    induced = [rng.random() < INDUCED_FRACTION for _ in tf_sites]
    dom_width = [rng.randrange(1_000, 3_000) for _ in domains]
    dom_strength = [math.exp(rng.gauss(0, 0.5)) for _ in domains]

    def bed(name, rows):
        with open(os.path.join(ref, name), "w") as out:
            for r in rows:
                out.write("\t".join(str(x) for x in r) + "\n")

    bed("blacklist.bed", [("chr22", c - 2_500, c + 2_500) for c in blacklist])
    bed("greenlist.bed", [("chr22", c - 500, c + 500, "Greenlist_region") for c in greenlist])
    bed("truth_hotspots.bed", [("chr22", c - 300, c + 300) for c in hotspots])
    bed("truth_tf_sites.bed", [("chr22", c - 100, c + 100, "induced" if i else "static", f"{s:.3f}")
                               for c, s, i in zip(tf_sites, tf_strength, induced)])
    bed("truth_domains.bed", [("chr22", c - w // 2, c + w // 2) for c, w in zip(domains, dom_width)])

    # GTF: transcripts starting at 40 % of TF sites, random strand
    with open(os.path.join(ref, "genes.gtf"), "w") as out:
        for i, c in enumerate(tf_sites):
            if rng.random() > 0.4:
                continue
            strand = rng.choice("+-")
            tss = c + 1
            s, e = (tss, min(L, tss + 5_000)) if strand == "+" else (max(1, tss - 5_000), tss)
            gid, tid = f"GENE{i:04d}", f"TX{i:04d}"
            attr = f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gid}";'
            for feat in ("gene", "transcript", "exon"):
                a_ = f'gene_id "{gid}"; gene_name "{gid}";' if feat == "gene" else attr
                out.write(f"chr22\tsim\t{feat}\t{s}\t{e}\t.\t{strand}\t.\t{a_}\n")

    rep_noise = {}
    for group, rep, run, kind in RUNS:
        key = (group, rep)
        if key not in rep_noise:
            rep_noise[key] = (
                [math.exp(rng.gauss(0, 0.2)) for _ in tf_sites],
                math.exp(rng.gauss(0, 0.1)),   # spike-in amount per reaction
            )
        noise, spike_noise = rep_noise[key]

        # target component
        target = None
        w_target = 0.0
        if kind.startswith("tf"):
            ws = [s * n * (INDUCTION if (kind == "tf_stim" and i) else 1.0)
                  for s, n, i in zip(tf_strength, noise, induced)]
            base = sum(s * n for s, n in zip(tf_strength, noise))
            w_target = W_TF * sum(ws) / base
            target = ("tf", Sampler(list(zip(tf_sites, ws))), 40)
        elif kind == "histone":
            w_target = W_HISTONE
            target = ("histone", Sampler(list(zip(zip(domains, dom_width), dom_strength))), None)

        comps = [
            ("background", W_BACKGROUND),
            ("greenlist", W_GREENLIST),
            ("hotspot", W_HOTSPOT * (3.0 if kind == "igg" else 1.0)),
            ("blacklist", W_BLACKLIST),
            ("mito", W_MITO),
            ("spike", W_SPIKE * spike_noise),
        ]
        if target:
            comps.append(("target", w_target))
        comp_sampler = Sampler(comps)

        r1p = os.path.join(fq, f"{group}_{run}_R1.fastq.gz")
        r2p = os.path.join(fq, f"{group}_{run}_R2.fastq.gz")
        counts = {c: 0 for c, _ in comps}
        with gzip.open(r1p, "wt", compresslevel=3) as o1, gzip.open(r2p, "wt", compresslevel=3) as o2:
            n = 0
            prev = None
            while n < a.pairs:
                if prev is not None and rng.random() < 0.08:
                    src, fs, fl, comp = prev       # PCR duplicate
                else:
                    comp = comp_sampler.pick(rng)
                    src = chrom
                    if comp == "background":
                        fl = frag_len(rng, "bg")
                        fs = rng.randrange(0, L - fl)
                    elif comp == "greenlist":
                        fl = frag_len(rng, "bg")
                        fs = int(rng.choice(greenlist) + rng.uniform(-450, 450) - fl / 2)
                    elif comp == "hotspot":
                        fl = frag_len(rng, "bg")
                        fs = int(rng.choice(hotspots) + rng.gauss(0, 150) - fl / 2)
                    elif comp == "blacklist":
                        fl = frag_len(rng, "bg")
                        fs = int(rng.choice(blacklist) + rng.gauss(0, 300) - fl / 2)
                    elif comp == "mito":
                        src, fl = chrm, frag_len(rng, "bg")
                        fs = rng.randrange(0, len(chrm) - fl)
                    elif comp == "spike":
                        src, fl = spike, frag_len(rng, "bg")
                        fs = rng.randrange(0, len(spike) - fl)
                    elif target[0] == "tf":
                        fl = frag_len(rng, "tf")
                        fs = int(target[1].pick(rng) + rng.gauss(0, target[2]) - fl / 2)
                    else:
                        c, w = target[1].pick(rng)
                        fl = frag_len(rng, "histone")
                        fs = int(c + rng.uniform(-w / 2, w / 2) - fl / 2)
                if fs < 0 or fs + fl > len(src):
                    continue
                frag = src[fs:fs + fl]
                if "N" in frag:
                    continue
                prev = (src, fs, fl, comp)
                counts[comp] += 1
                if rng.random() < 0.5:
                    frag = revcomp(frag)
                r1 = (frag + ADAPTER_R1 + "A" * READ_LEN)[:READ_LEN]
                r2 = (revcomp(frag) + ADAPTER_R2 + "A" * READ_LEN)[:READ_LEN]
                r1, r2 = mutate(r1, rng), mutate(r2, rng)
                q = "I" * READ_LEN
                name = f"@{group}_{run}_{n}"
                o1.write(f"{name}/1\n{r1}\n+\n{q}\n")
                o2.write(f"{name}/2\n{r2}\n+\n{q}\n")
                n += 1
        mix = ", ".join(f"{c} {100 * v / a.pairs:.1f}%" for c, v in counts.items())
        print(f"  {group} {run}: {a.pairs} pairs ({mix})")


if __name__ == "__main__":
    main()
