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
