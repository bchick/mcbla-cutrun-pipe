# Outputs

All paths are relative to the Snakemake working directory. Logs mirror this layout under `logs/`. `<lib>` is `<group>_R<replicate>`, and `<caller>` is `macs2` or `seacr`.

## Reference (`results/reference/`)

| file | content |
|---|---|
| `genome.fa`, `.fai`, `chrom.sizes` | host FASTA, linked or decompressed, plus its index |
| `host.chrom.sizes` | contigs without mito or spike-in; used for background windows and bedGraphs |
| `bowtie2/`, `bowtie2_combined/` | indexes built by the pipeline when none is given (host, or host + spike-in) |
| `combined.fa` | host + spike-in FASTA, with spike-in contigs renamed `<contig_prefix><name>` |
| `tss.bed` | TSSs from the GTF, used for the TSS heatmaps |

## BAMs (`results/bam/`)

| file | content |
|---|---|
| `<lib>.final.bam` (+ `.bai`) | filtered (MAPQ, proper pairs, mito, spike-in, blacklist) and duplicate-marked. Duplicates are **removed** for controls when `duplicates.remove_control` is set, and for targets when `remove_target` is set. |
| `sized/<lib>.bam` | fragments ≤ `fragments.peak_max_size`, only when it is > 0 |
| `merged/<group>.bam`, `merged_sized/<group>.bam` | replicates merged, for group peaks and as MACS2 controls |
| `igg_pool/pool.bam` | all control libraries merged, for hotspots, the fingerprint JS distance and the gate for groups without a control |

## QC (`results/qc/`)

| file | content |
|---|---|
| `alignment_qc_report.tsv` | fastq mode: raw, trimmed, aligned, mito, spike-in, blacklist, final, duplicate %, mean fragment size, TF fraction, NRF, PBC1, PBC2 |
| `qc_summary.tsv` | per library: PASS/WARN/FAIL per metric (`qc.thresholds`), an overall flag and the metrics not passing |
| `igg_enrichment.tsv` | per library and group peak set: median fold over IgG, fraction of peaks ≥ `igg.min_fold`, fraction on IgG hotspots, peaks kept |
| `spikein_summary.tsv` | spike-in fragments and % of host fragments (when there is a spike-in) |
| `flagstat/`, `markdup/`, `filter_stats/`, `complexity/`, `fragment_sizes/`, `frip/`, `spikein/`, `fastqc/` | per-library inputs to the tables above |
| `deeptools/` | fragment-size plot and table; Spearman correlation and PCA (500 bp bins, blacklist excluded); fingerprint and its metrics (JS distance to the pooled IgG) |
| `multiqc/multiqc_report.html` | everything above in one report |

## Peaks (`results/peaks/`)

| file | content |
|---|---|
| `<caller>/individual/<lib>.raw.bed` | per-library peaks after the blacklist (narrowPeak columns; broad and SEACR peaks have summit −1 or the centre of the maximum block) |
| `<caller>/merged/<group>.raw.bed` | peaks called on the merged-replicate BAM |
| `<caller>/{individual,merged}/<name>.bed` | **the peaks used downstream**: `.raw.bed` after the IgG gate (identical to it unless `igg.hotspot_filter` is on) |
| `<caller>/.../<name>.igg.tsv`, `.igg_summary.tsv` | per-peak target/IgG fragments, fold, hotspot overlap, pass and kept flags (with IgG libraries) |
| `macs2/.../<name>_peaks.xls` and native MACS2 files | the MACS2 output as written |
| `<caller>/{individual,merged}/peak_summary.tsv` | peaks, FRiP (fragments, duplicates excluded), median width |
| `consensus/groups/<group>.bed` | peaks in ≥ min(`peaks.min_overlap`, n replicates) replicates of the group |
| `consensus/<target>.consensus.bed` | union of the target's group consensus sets (`peaks: consensus`) |
| `consensus/<target>.merged.bed` | union of the target's merged-BAM peaks (`peaks: merged`) |
| `consensus/<target>.idr.bed` | union of IDR peaks, with merged peaks for single-replicate and broad groups (`peaks: idr`) |
| `consensus/consensus_summary.tsv` | peak counts and median widths of the sets above |
| `idr/` | relaxed peaks, IDR per replicate pair, `<group>_idr.narrowPeak` and `idr_summary.tsv` (`modules.idr`) |
| `results/igg/hotspots.bed`, `background_windows.bed` | pooled-IgG hotspots (plus `igg.hotspot_bed`) and the random windows used for depth-free scaling |

## Signal (`results/bigwig/`, `results/normalization/`)

| file | content |
|---|---|
| `bigwig/<method>/<lib>.bw` | per-library track: `depth` (CPM/RPGC), `greenlist` or `spikein` (CPM-equivalent) |
| `bigwig/<method>/groups/<group>.bw` | replicate average (bigwigAverage) |
| `normalization/<method>/size_factors.tsv` | SampleID, group, target, norm_group, size_factor, reads, scale_factor (the bamCoverage `--scaleFactor`) |
| `normalization/greenlist/greenlist_counts.tsv` | fragments per greenlist region per library |

## Analyses

| file | content |
|---|---|
| `diff/<target>/dba_counted.rds` | counted DiffBind object shared by every normalization |
| `diff/<target>/<method>/summary.tsv` | gained, lost and significant counts per contrast |
| `diff/<target>/<method>/tables/<label>_{all,sig}.tsv` | full and FDR-significant DiffBind reports |
| `diff/<target>/<method>/size_factors.tsv`, `plots.pdf`, `dba_analyzed.rds` | normalization factors used; correlation, PCA, MA and bar plots; analysed object |
| `diff/contrasts/<label>_{gained,lost}.bed` | significant peaks of the `diff.normalization` run |
| `diff/<target>/normcheck/norm_comparison.tsv`, `norm_verdict.tsv`, `norm_comparison_barplot.pdf` | normalization sensitivity per contrast and method, with a recommendation |
| `heatmaps/<target>/{peaks,tss}_heatmap.png` (+ `_matrix.gz`) | group tracks (and their IgG) over the target's peaks and over TSSs |
| `heatmaps/contrasts/<label>_heatmap.png` | the contrast's two groups over its gained and lost peaks |
| `annotate/<set>.annotation.tsv`, `annotation_summary.tsv`, `annotation_plots.pdf` | ChIPseeker annotation per target peak set and per contrast direction |
| `motifs/<target>/`, `motifs/contrasts/<label>_{gained,lost}/` | HOMER `knownResults.txt` and `knownResults.html` (plus `homerResults/` with `denovo: true`) |
