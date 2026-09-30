# Test dataset

A synthetic paired-end CUT&RUN dataset on a 4 Mb window of chr22, built by
`scripts/build_test_data.sh` (standard-library Python, ~25 MB, gitignored).
See the docstring of `scripts/simulate_cutrun.py` for the model.

| group | libraries | notes |
|---|---|---|
| tfA_ctrl | 2 | TF, IgG control |
| tfA_stim | 3 (rep 3 in two runs) | TF with 30 % of sites induced 4x: a global increase in binding |
| k27_ctrl | 1 | histone, broad, **no control** (exercises the IgG gate on the pooled IgG) |
| IgG | 2 | background, greenlist regions and 15 sticky hotspots |

Every library also carries greenlist regions at a constant absolute level and
spike-in (E. coli-like) carry-over at a constant amount per reaction.

`config/config.yaml` switches every option on (spike-in alignment, IgG
control + gate, SEACR, IDR, sub-nucleosomal peak calling, all three
normalizations, diff with greenlist, normcheck, heatmaps, annotation, HOMER).

    pixi run build-test   # build the data
    pixi run test         # align -> peaks -> consensus (test_core)
    pixi run test-all     # everything, then scripts/check_results.py

`check_results.py` checks the known answers: peaks on the simulated TF sites,
the IgG gate removing hotspot peaks, greenlist and spike-in size factors
agreeing on the global shift, depth normalization's false "lost" peaks
removed by greenlist, and normcheck flagging the contrast.

## Head-to-head with nf-core/cutandrun

`nfcore/` runs nf-core/cutandrun 3.2.2 (MACS2, IgG control, spike-in) on the
same simulated reads and scores both pipelines against the simulated truth.
It needs Nextflow and Docker (`NFCORE_PROFILE=singularity` to switch) and the
`test-all` outputs:

    pixi run test-all
    pixi run test-nfcore      # nf-core run (~20 min at 8 cores, then -resume) + compare
    pixi run compare-nfcore   # comparison only, on existing outputs

`nfcore/compare_nfcore.py` writes `results/nfcore_compare/{metrics.tsv,report.md}`.
Only these checks can fail it (tolerances in the script docstring):

- precision and F1 of the tfA consensus vs `truth_tf_sites.bed`, compared
  with nf-core's pooled tfA_ctrl + tfA_stim consensus;
- k27 domain recall, and k27 peaks on IgG hotspots only (must be no more
  than nf-core's);
- per-library duplicate rate vs the simulated 8%.

Peak counts, the peak-set Jaccard, FRiP and the spike-in stim/ctrl ratio are
reported only. Two things the runner has to work around in 3.2.2: it cannot
leave one target uncontrolled while others use IgG, so k27_ctrl is called
against IgG there (ours calls it without a control and applies the IgG gate);
and it runs under Nextflow 24.04.4 (`NXF_VER`), because its Trim Galore module
does not parse under Nextflow 25.
