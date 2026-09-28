# mcbla-cutrun-pipe

A Snakemake workflow for paired-end **CUT&RUN**. It takes FASTQ files (or filtered BAMs) through to QC, peaks, normalized signal tracks and differential binding.

- **IgG controls are optional.** You can use them in three independent ways: as the peak-calling control, as an enrichment filter, and for QC.
- **Three normalizations.** Depth is the default. Greenlist and spike-in are available for experiments that change total binding.
- **Downstream analyses are opt-in.** DiffBind, a normalization-sensitivity check, heatmaps, peak annotation and HOMER motifs each switch on separately.

The processing defaults are those of the McBla lab CUT&RUN scripts in `mcf7_project`, which lab benchmarks chose (peak callers, aligners, normalization). The infrastructure follows [mcbla-bulkatac-pipe](https://github.com/bchick/mcbla-bulkatac-pipe): pixi launcher, pinned per-step conda envs, Slurm and local profiles, a validated samplesheet and config, a synthetic test dataset, and CI.

```
FASTQ ─ cutadapt ─ bowtie2 (host [+ spike-in]) ─ MAPQ/pair filter ─ mito/blacklist ─ mark duplicates
                                                         │
             ┌───────────────────────────────────────────┼────────────────────────────┐
         QC (FastQC, fragment sizes,             peaks: MACS2 [-c IgG] | SEACR     signal: depth | greenlist
         TF fraction, FRiP, correlation,         → blacklist → IgG gate            | spike-in bigWigs
         fingerprint vs IgG, qc_summary,         → group consensus (≥ k reps)      (per library, per group)
         MultiQC)                                → per-target peak sets [IDR]
                                                         │
                                  opt-in: DiffBind (per target) · normcheck · heatmaps · ChIPseeker · HOMER
```

## Quick start

**On the Salk lab server**, `pixi run init` asks which genome you are using, whether you have IgG controls, which spike-in and which blacklist. It then writes a project config that points at the shared references listed in `/data/resource/manifest.yaml`:

```bash
pixi install
pixi run init                      # interactive; or pass --dir --genome --igg --spikein --blacklist
pixi run init --list               # the genomes, blacklists and spike-ins on offer
```

It prints the dry-run command for the new project. Agents running the pipeline for someone should follow [AGENTS.md](AGENTS.md).

**Anywhere else:**

```bash
git clone https://github.com/bchick/mcbla-cutrun-pipe.git
cd mcbla-cutrun-pipe
pixi install                       # Snakemake + plugins only; tools come from workflow/envs/

# 1. describe your libraries: config/samples.tsv (and config/contrasts.tsv for diff)
# 2. point config/config.yaml (or a copy of config/salk_example.yaml) at your references
pixi run snakemake -n                                           # dry run
pixi run snakemake --profile profiles/slurm                     # cluster
pixi run snakemake --profile profiles/local --cores 32          # one machine
```

To keep a project outside the repo, give its directory and a config that holds only the keys you change. The repo's `config/config.yaml` supplies the rest.

```bash
pixi run snakemake -s workflow/Snakefile --directory /path/to/project \
    --configfile /path/to/project/project.yaml --profile profiles/slurm
```

Try it first on the synthetic test dataset, which takes a few minutes:

```bash
pixi run build-test     # ~25 MB of simulated reads on a 4 Mb chr22 window
pixi run test           # align -> peaks -> consensus
pixi run test-all       # every module, then .test/scripts/check_results.py
pixi run test-init      # `pixi run init` against a fixture manifest, every choice dry-run
```

The config, samplesheet and contrasts are checked before any job runs. Mistakes stop the run with a plain-language message, for example `contrast X: FOS_HRG60 (FOS) and K27ac_US (H3K27ac) are different targets`.

## Samplesheet

TSV or CSV. The first five columns follow **nf-core/cutandrun**.

| column | required | meaning |
|---|---|---|
| `group` | yes | experimental group, antibody × condition (e.g. `FOS_HRG60`). The unit of merging, consensus and contrasts. |
| `replicate` | yes | biological replicate (1, 2, …). Libraries are named `<group>_R<replicate>`. |
| `fastq_1`, `fastq_2` | fastq mode | paired-end reads. Rows sharing group and replicate are sequencing runs of one library and are concatenated. |
| `control` | no | the group holding this group's IgG libraries. Leave it empty for the IgG rows themselves. |
| `target` | no | antibody (`FOS`, `H3K27ac`, `IgG`). Groups of one target share a consensus peak set, one DiffBind object and, by default, greenlist size factors. Defaults to the group name. |
| `target_type` | no | `tf` (default) or `histone`. Only `tf` libraries are flagged on the sub-nucleosomal fragment fraction and size-selected for peak calling. |
| `peak_mode` | no | `narrow` (default) or `broad` (MACS2 `--broad`, e.g. H3K27me3, H3K36me3). |
| `batch` | no | prep or sequencing batch. Used for the IgG batch-mismatch warning, `greenlist_group: target_batch` and the `batch: true` covariate. |
| `bam` | bam mode | filtered, coordinate-sorted BAM. |
| `spikein_reads` | no | bam mode: spike-in fragments per library, which enables spike-in normalization. |

```
group      replicate  fastq_1                 fastq_2                 control  target   target_type
FOS_US     1          FOS_US_r1_R1.fastq.gz   FOS_US_r1_R2.fastq.gz   IgG      FOS      tf
FOS_US     2          ...                                             IgG      FOS      tf
FOS_HRG60  1          ...                                             IgG      FOS      tf
IgG        1          IgG_r1_R1.fastq.gz      IgG_r1_R2.fastq.gz               IgG
```

Contrasts (`group1  group2  label`, with log2FC = group1 / group2) compare two groups of the **same target**. Each group in a contrast needs at least 2 replicates.

## IgG controls

Every IgG option is independent and can be switched off. A group is treated as a control if it is named in another group's `control` column or if its `target` is `IgG`. Control groups get BAMs, bigWigs and QC, but no peaks or contrasts. By default, duplicates are removed from controls and only marked on targets (`duplicates:`).

| option | default | what it does |
|---|---|---|
| `igg.as_control` | on | MACS2 `-c <merged IgG of the group>`; SEACR in control mode. Groups without a `control` are called without one. |
| `igg.qc` | on | For every peak set, reports the **fold over IgG** (median and fraction ≥ `min_fold`) and the overlap with pooled-IgG hotspots. The results go to `qc_summary.tsv` and MultiQC. |
| `igg.hotspot_filter` | off | Drops peaks whose fold over IgG is below `igg.min_fold` (2×). The unfiltered peaks stay in `<name>.raw.bed`. With `remove_hotspot_overlaps`, it also drops peaks overlapping pooled-IgG MACS2 peaks or `hotspot_bed`. |

**How the fold is computed.** Target and IgG fragment centres are counted in each peak and in random background windows away from peaks. Each library is divided by its own background density, so the fold does not depend on library depth or FRiP. This is the design of the lab's IgG reality gate (mcf7 analysis 39). The reference IgG is the group's `control`, or all IgG libraries pooled when the group has none.

**When to use which option:**

- **IgG from the same batch:** keep `as_control` on.
- **IgG only from another batch:** the pipeline warns. In the mcf7 data, a mismatched IgG *added* artefactual peaks, and a pooled IgG made no difference to peak calls. Consider `as_control: false` with `hotspot_filter: true`: the IgG then acts only as a filter for sticky regions.
- **No IgG at all:** run everything without it. You can give an external hotspot BED as `igg.hotspot_bed`.

## Normalization

`normalization.bigwig_methods` picks which tracks are written, and `diff.normalization` picks what DiffBind uses.

| method | how | use for |
|---|---|---|
| `depth` (default) | bamCoverage CPM (or RPGC); DiffBind library size | replicate QC and PCA, and contrasts where total binding does not change |
| `greenlist` | DESeq2 median-of-ratios on the 868-region CUT&RUN greenlist (de Mello et al. 2024; hg38 bundled, `reference.greenlist`) | contrasts with a **global shift**: stimulated vs unstimulated, degrader vs vehicle |
| `spikein` | spike-in fragments (E. coli carry-over or added yeast/fly DNA) | the same, when the spike-in is reliable |

- **Size factors are estimated within a target** (`greenlist_group: target_batch` estimates them within target × batch). Factors estimated over more samples than are being compared can differ substantially. In the mcf7 SD51 data, project-wide and subset factors correlated at only r = 0.85.
- **Greenlist and spike-in tracks share CPM units.** They are scaled to be *CPM-equivalent*: identical to the CPM track when there is no global shift. They can therefore be compared by eye with the depth tracks. Every track is written under `results/bigwig/<method>/` so methods are not mixed by accident.
- **normcheck** (`normcheck.run`) re-runs every contrast under the other methods on the same counts. It flags a contrast when its gained or lost counts move by more than 20 % (and by at least 10 peaks), or when the direction of the net change flips. This is the decision rule of mcf7 analyses 20 and 22, where TMM flipped the sign of the MEKi contrast and greenlist was chosen for global-shift contrasts.

## Peaks

- **MACS2** is the default caller (`-f BAMPE --keep-dup all -q 0.05`, or `--broad` per group) and won the lab benchmark (mcf7 analysis 19).
- **SEACR** is opt-in (`peaks.seacr.run`, stringent by default). It works well with an IgG control. Without one it applies a numeric threshold, which gave degenerate peak sets in the benchmark, and the pipeline warns about this.
- **Blacklist filtering** applies to every peak set, and **FRiP counts fragments** with duplicates excluded.
- **Consensus per group** keeps merged peak regions that are supported by at least `min(peaks.min_overlap, n replicates)` replicates, retaining the full peak extents. Groups come from the samplesheet, so treatment arms at the same time point are never pooled.
- **Per-target peak sets** feed the analyses:
  - `consensus`: the union of the group consensus sets;
  - `merged`: peaks called on merged-replicate BAMs;
  - `idr`: IDR-reproducible peaks (with `modules.idr`).
- **Sub-nucleosomal calling (optional):** `fragments.peak_max_size: 120` calls TF peaks on fragments of 120 bp or less only, and size-selects their IgG the same way.

## Analyses (opt-in)

Each analysis runs only with `<section>.run: true`, and you then set its `peaks:` to `consensus`, `merged`, `idr`, `individual` (diff only) or a BED path. Outputs from earlier steps are reused, so switching an analysis on later only runs that analysis.

| section | what |
|---|---|
| `diff` | DiffBind (DESeq2 backend), one object per target. Design `~group` (or `~batch + group`). Produces tables, gained/lost BEDs, MA, PCA and correlation plots. Normalization: `depth`, `greenlist` or `spikein`. |
| `normcheck` | the same contrasts under `normcheck.methods` (`depth`, `greenlist`, `spikein`, `csaw` background bins), plus a verdict per contrast |
| `heatmaps` | deepTools heatmaps of replicate-averaged tracks (target groups plus their IgG) over peaks, TSSs and each contrast's gained/lost peaks |
| `annotate` | ChIPseeker, using a TxDb built from `reference.gtf` (works for any genome) |
| `motifs` | HOMER known motifs (and de novo with `denovo: true`). Target peaks are tested against a GC-matched genomic background, and gained/lost peaks against the target's peak set. |

## Input modes

- `input_mode: fastq` runs the full pipeline.
- `input_mode: bam` starts from filtered BAMs given in the samplesheet `bam` column, for example the `*_markdup.bam` files of an earlier run. Controls are deduplicated when `duplicates.remove_control` is on. Spike-in normalization then comes from the `spikein_reads` column.

## Profiles and environments

- `profiles/slurm`: set `slurm_partition` and `slurm_account`.
- `profiles/local`: a single machine.
- **Tools** are pinned per step in `workflow/envs/*.yaml` (bowtie2 2.5.2, cutadapt 4.6, samtools 1.13, MACS2 2.2.9.1, deepTools 3.5.4, DiffBind 3.10). These are the versions used for the mcf7 analyses.
- **Shared env store:** `SNAKEMAKE_CONDA_PREFIX` sets where envs are built. Point it at a shared lab directory so each lab member does not rebuild them.

## Outputs

See [docs/outputs.md](docs/outputs.md). The main files:

- `results/qc/multiqc/multiqc_report.html`
- `results/qc/qc_summary.tsv`: PASS/WARN/FAIL per library for fragments, alignment, mito, duplicates, FRiP, TF fraction and IgG enrichment.
  The flags are reported, never enforced. Expect FRiP to FAIL for unstimulated samples of an inducible TF (e.g. FOS at baseline has few peaks); judge those against their stimulated groups.
- `results/peaks/consensus/<target>.{consensus,merged,idr}.bed`
- `results/bigwig/<method>/groups/<group>.bw`
- `results/diff/<target>/<method>/`

## Defaults and their evidence

See [docs/defaults_rationale.md](docs/defaults_rationale.md).

## Known gaps

- Single-end reads are not supported.
- There are no IDR pseudo-replicates, so single-replicate groups get no reproducibility estimate.
- The greenlist is bundled for hg38 only; give a BED for other genomes.
- Spike-in normalization from an external count table works in bam mode only.

## Citation

If you use this workflow, cite it (see `CITATION.cff`) and the tools it runs:

- Bowtie2 (Langmead & Salzberg 2012)
- cutadapt (Martin 2011)
- SAMtools (Danecek et al. 2021)
- BEDTools (Quinlan & Hall 2010)
- MACS2 (Zhang et al. 2008)
- SEACR (Meers et al. 2019)
- deepTools (Ramírez et al. 2016)
- IDR (Li et al. 2011)
- DiffBind (Ross-Innes et al. 2012; Stark & Brown)
- DESeq2 (Love et al. 2014)
- csaw (Lun & Smyth 2016)
- ChIPseeker (Yu et al. 2015)
- HOMER (Heinz et al. 2010)
- MultiQC (Ewels et al. 2016)
- the CUT&RUN greenlist (de Mello et al. 2024)
