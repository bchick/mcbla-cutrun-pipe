# Why the defaults are what they are

Most defaults come from the McBla lab's CUT&RUN processing (`mcf7_project/scripts/alignment/1.1_cutrun_*`, `scripts/cutrun/1.2–1.8`) and from the lab benchmarks in `mcf7_project/analyses/19–22`. Where this pipeline departs from those scripts, the reason is given here.

## Alignment and filtering

- **bowtie2 `--very-sensitive-local -I 10 -X 700 --dovetail --no-mixed --no-discordant`.** These are the lab settings, and benchmark 21 found no reason to switch aligners.
  - `--dovetail` keeps pairs whose mates extend past each other. This is common for CUT&RUN fragments shorter than the read length.
  - `-I 10` admits short TF footprints.
- **MAPQ ≥ 20, proper pairs, `-F 2828`, mito removed, ENCODE blacklist.** These are the lab settings.
- **Duplicates are marked on targets, not removed.** In CUT&RUN, MNase cuts at a bound site recur at the same position, so real signal looks like PCR duplicates. The lab decision was that deduplication over-removes signal. MACS2 runs with `--keep-dup all`.
- **Duplicates are removed from IgG controls.** IgG libraries are low-complexity, and their duplicates are artefacts. nf-core/cutandrun does the same by default.
- **FRiP counts fragments, duplicates excluded.** The lab script counted reads with duplicates included, which inflates FRiP for low-complexity libraries.

## Peaks

- **MACS2 (`-f BAMPE --keep-dup all -q 0.05`).** Benchmark 19 compared peak callers on the lab data and kept MACS2.
- **SEACR is opt-in.** SEACR's numeric mode (top 1 % of signal blocks, used without IgG) gave degenerate peak sets in benchmark 19. It is offered because it is the reference caller when a matched IgG exists (Meers et al. 2019).
- **`-g hs`.** This is the lab setting. The deepTools effective genome size (2 913 022 398) is used only for RPGC tracks.
- **Consensus requires ≥ 2 replicates and keeps full peak extents.**
  - The lab `1.5` script intersected replicate pairs. It also grouped by the `COND_TIME` part of the filename, which pooled treatment arms: SD46 `veh` and `meki` at HRG_60m were merged.
  - Here, groups come from the samplesheet, and a region needs peaks from k distinct replicates.

## IgG

- **The IgG is a control when given, but nothing requires it.** The mcf7 experience was that a pooled IgG made no difference to peak calls and that a mismatched IgG added artefacts. The pipeline therefore:
  - uses an IgG only when the samplesheet names one;
  - warns when it comes from another batch;
  - offers the IgG as a filter instead of a control.
- **The enrichment gate is depth-free (analysis 39).** Raw target/IgG count ratios are really ratios of library depth. Dividing each library by its own density over random background windows removes depth without needing a size factor, and it is not biased by the target's FRiP the way CPM is.
  - Target and IgG are counted with the same instrument (fragment centres, same MAPQ).
  - The 2× threshold is the analysis 39 value.
- **Hotspot overlaps are reported, not removed, by default.** Real peaks at open promoters often overlap IgG signal. The fold gate is the discriminating test, and hotspot removal is available with `remove_hotspot_overlaps`.

## Normalization

- **Depth is the default; greenlist is the recommendation for global shifts.** Analyses 20 and 22 compared normalizations on the lab's AP-1, BAF and H3K27ac data:
  - library-size (depth) normalization is fine for replicate QC and for within-condition comparisons;
  - for global-shift contrasts (stimulation vs unstimulated, degrader vs vehicle), greenlist normalization was chosen;
  - TMM flipped the sign of the MEKi contrast;
  - spike-in (E. coli carry-over) was unstable in the lab data.
- **normcheck makes the decision explicit for every contrast.** It flags a contrast when its gained or lost counts change by more than 20 % and by at least 10 peaks, or when the direction of the net change flips.
- **Greenlist size factors are computed within a target, or within target × batch.** For SD51, project-wide and subset size factors correlated at r = 0.85. Factors should be computed on the libraries being compared.
- **Greenlist and spike-in tracks are CPM-equivalent.** The scale is 1e6 / (size factor × geometric-mean mapped reads of the normalization group). Without a global shift this reduces to CPM, so all methods share one unit.
  - The lab `bw_greenlist/` tracks were raw coverage divided by the size factor, and could not be compared with `bw/`.
  - The default depth track is therefore CPM rather than RPGC, so its units match the other methods.
- **Greenlist size factors use DESeq2 median-of-ratios.** They are computed in Python with the same formula as `DESeq2::estimateSizeFactorsForMatrix`: rows with a zero in any library are dropped, and each library's factor is the median of its ratio to the row geometric mean.

## Differential binding

- **DiffBind 3.10 with the DESeq2 backend and an explicit design.** These are the mcf7 versions. Contrasts are added under a design (`~Condition`, or `~Factor + Condition` with `batch: true`).
  - If contrasts are given as group masks instead, DiffBind 3.x falls back to a legacy per-contrast normalization that ignores `dba.normalize()`. Every normalization then gives the same p-values; mcbla-bulkatac-pipe found and fixed this.
  - `check_size_factors()` stops the run if DESeq2 did not use the stored factors.
- **Greenlist and spike-in factors enter through `dba.normalize(library = sf, normalize = DBA_NORM_LIB)`,** as in the lab's `apply_greenlist_to_dba()`.
- **Each contrasted group needs ≥ 2 replicates.** DiffBind refuses otherwise, so the pipeline checks this before any job runs.

## Tool versions

All tools are pinned per step to the versions used for the mcf7 analyses (`pipeline-tools/pixi.toml`): bowtie2 2.5.2, cutadapt 4.6, samtools 1.13, bedtools 2.30.0, MACS2 2.2.9.1, deepTools 3.5.4, IDR 2.0.4.2, MultiQC 1.17, DiffBind 3.10 / DESeq2 1.40 (Bioconductor 3.17).

The MACS2 build is pinned (`py311hdad781d_1`) because newer bioconda rebuilds fail at import on current glibc. These pinned versions do not co-solve in one environment, so each step has its own.
