# Shared helpers: config and samplesheet loading, validation, control (IgG)
# resolution, peak-set and normalization resolution, and target lists.
#
# All Python functions live here, so the rule files contain only rules.

import os
import re
import sys
import glob
import itertools
from pathlib import Path

import pandas as pd
from snakemake.utils import validate

# ---------------------------------------------------------------------------
# Config + samplesheet
# ---------------------------------------------------------------------------
validate(config, schema="../schemas/config.schema.yaml")

INPUT_MODE = config["input_mode"]
MODULES = config["modules"]
IGG = config["igg"]
PK = config["peaks"]
NORM = config["normalization"]


def _warn(msg):
    print(f"[mcbla-cutrun-pipe] WARNING: {msg}", file=sys.stderr)


def _read_table(path):
    return pd.read_csv(path, sep=None, engine="python", dtype=str, comment="#").fillna(
        ""
    )


OPTIONAL_COLUMNS = (
    "fastq_1",
    "fastq_2",
    "control",
    "target",
    "target_type",
    "peak_mode",
    "batch",
    "bam",
    "spikein_reads",
)

samples_raw = _read_table(config["samples"])
for col in OPTIONAL_COLUMNS:
    if col not in samples_raw.columns:
        samples_raw[col] = ""
samples_raw["target"] = [
    t if t else g for t, g in zip(samples_raw["target"], samples_raw["group"])
]
samples_raw["target_type"] = [t or "tf" for t in samples_raw["target_type"]]
samples_raw["peak_mode"] = [m or "narrow" for m in samples_raw["peak_mode"]]
validate(samples_raw, schema="../schemas/samples.schema.yaml")

if INPUT_MODE == "fastq":
    _nofq = samples_raw[(samples_raw["fastq_1"] == "") | (samples_raw["fastq_2"] == "")]
    if len(_nofq):
        raise ValueError(
            "input_mode=fastq needs paired-end fastq_1 and fastq_2 for every row; "
            f"missing for group(s): {', '.join(sorted(set(_nofq['group'])))}"
        )
else:
    _nobam = samples_raw[samples_raw["bam"] == ""]
    if len(_nobam):
        raise ValueError(
            "input_mode=bam needs a `bam` column value for every row; missing "
            f"for group(s): {', '.join(sorted(set(_nobam['group'])))}"
        )

# nf-core/cutandrun names each library <group>_R<replicate>; rows that share
# group + replicate are sequencing runs of one library (merged before trimming).
samples_raw["lib"] = [
    f"{g}_R{r}" for g, r in zip(samples_raw["group"], samples_raw["replicate"])
]
if INPUT_MODE == "bam":
    _multi = samples_raw["lib"][samples_raw["lib"].duplicated()]
    if len(_multi):
        raise ValueError(
            "input_mode=bam takes one BAM per library; repeated group and "
            f"replicate: {', '.join(sorted(set(_multi)))}"
        )

LIB_COLUMNS = ("group", "replicate", "control", "target", "target_type", "peak_mode")
_per_lib = samples_raw.groupby("lib", sort=False).agg(
    {
        **{c: lambda x: ";".join(sorted(set(x))) for c in LIB_COLUMNS},
        "batch": lambda x: ";".join(sorted(set(x))),
        "bam": "first",
        "spikein_reads": "first",
    }
)
for col in LIB_COLUMNS + ("batch",):
    bad = _per_lib[_per_lib[col].str.contains(";")]
    if len(bad):
        raise ValueError(
            f"Samplesheet: runs of the same library disagree on '{col}': "
            f"{', '.join(bad.index)}"
        )
LIBRARIES = _per_lib
LIBS = list(LIBRARIES.index)
GROUPS = list(dict.fromkeys(LIBRARIES["group"]))
LIBS_BY_GROUP = {g: list(LIBRARIES.index[LIBRARIES["group"] == g]) for g in GROUPS}

# Group-level attributes must agree across the group's libraries.
GROUP_INFO = {}
for _g in GROUPS:
    _rows = LIBRARIES.loc[LIBS_BY_GROUP[_g]]
    for col in ("control", "target", "target_type", "peak_mode"):
        if _rows[col].nunique() > 1:
            raise ValueError(
                f"Samplesheet: libraries of group '{_g}' disagree on '{col}' "
                f"({', '.join(sorted(set(_rows[col])))}); it must be the same for "
                "every replicate of a group."
            )
    GROUP_INFO[_g] = _rows.iloc[0][["control", "target", "target_type", "peak_mode"]]

# ---------------------------------------------------------------------------
# Controls (IgG)
# ---------------------------------------------------------------------------
# A control group is any group named in another group's `control` column, or
# whose target is IgG. Control groups get BAMs, bigWigs and QC but no peaks,
# consensus sets or contrasts.
_named_controls = {GROUP_INFO[g]["control"] for g in GROUPS} - {""}
_unknown = _named_controls - set(GROUPS)
if _unknown:
    raise ValueError(
        "Samplesheet `control` names group(s) that are not in the samplesheet: "
        + ", ".join(sorted(_unknown))
    )
CONTROL_GROUPS = [
    g
    for g in GROUPS
    if g in _named_controls or GROUP_INFO[g]["target"].lower() == "igg"
]
for _g in CONTROL_GROUPS:
    if GROUP_INFO[_g]["control"]:
        raise ValueError(
            f"Samplesheet: '{_g}' is a control group (IgG, or named as another "
            "group's control) but has a control of its own; leave its `control` empty."
        )
TARGET_GROUPS = [g for g in GROUPS if g not in CONTROL_GROUPS]
if not TARGET_GROUPS:
    raise ValueError("Samplesheet: every group is a control group; nothing to call.")
CONTROL_LIBS = [l for g in CONTROL_GROUPS for l in LIBS_BY_GROUP[g]]
TARGET_LIBS = [l for g in TARGET_GROUPS for l in LIBS_BY_GROUP[g]]
HAS_IGG = bool(CONTROL_LIBS)

TARGETS = list(dict.fromkeys(GROUP_INFO[g]["target"] for g in TARGET_GROUPS))
GROUPS_BY_TARGET = {
    t: [g for g in TARGET_GROUPS if GROUP_INFO[g]["target"] == t] for t in TARGETS
}


def lib_group(lib):
    return LIBRARIES.loc[lib, "group"]


def is_control_lib(lib):
    return lib in CONTROL_LIBS


def group_control(group):
    """The control group used for peak calling (igg.as_control), or None."""
    c = GROUP_INFO[group]["control"]
    return c if c and IGG["as_control"] else None


def _batches(libs):
    return {LIBRARIES.loc[l, "batch"] for l in libs} - {""}


for _g in TARGET_GROUPS:
    _c = group_control(_g)
    if not _c:
        continue
    _tb, _cb = _batches(LIBS_BY_GROUP[_g]), _batches(LIBS_BY_GROUP[_c])
    if _tb and _cb and not (_tb & _cb):
        _warn(
            f"group '{_g}' (batch {','.join(sorted(_tb))}) uses control '{_c}' "
            f"from another batch ({','.join(sorted(_cb))}). A mismatched IgG can "
            "add artefactual peaks; prefer a same-batch IgG or igg.as_control: false."
        )
if IGG["as_control"] and HAS_IGG:
    _noctl = [g for g in TARGET_GROUPS if not GROUP_INFO[g]["control"]]
    if _noctl:
        _warn(
            "igg.as_control is true, but these groups have no `control` and are "
            "called without one: " + ", ".join(_noctl)
        )

# IgG-based peak QC / filtering needs IgG libraries, or (filter only) a BED.
IGG_QC = bool(IGG["qc"]) and HAS_IGG
IGG_FILTER = bool(IGG["hotspot_filter"])
if IGG_FILTER and not HAS_IGG and not IGG.get("hotspot_bed"):
    raise ValueError(
        "igg.hotspot_filter is true, but the samplesheet has no IgG/control "
        "libraries and igg.hotspot_bed is empty. Add IgG libraries, give a "
        "hotspot BED, or set igg.hotspot_filter: false."
    )
IGG_FOLD = HAS_IGG and (IGG_QC or IGG_FILTER)
USE_HOTSPOTS = HAS_IGG or bool(IGG.get("hotspot_bed"))

# ---------------------------------------------------------------------------
# Contrasts
# ---------------------------------------------------------------------------
contrasts = _read_table(config["contrasts"]) if config.get("contrasts") else None
if contrasts is not None and len(contrasts):
    validate(contrasts, schema="../schemas/contrasts.schema.yaml")
    _unknown = (set(contrasts["group1"]) | set(contrasts["group2"])) - set(GROUPS)
    if _unknown:
        raise ValueError(f"contrasts.tsv references unknown groups: {sorted(_unknown)}")
    CONTRASTS = list(contrasts["label"])
    if len(set(CONTRASTS)) != len(CONTRASTS):
        raise ValueError("contrasts.tsv: labels must be unique")
    CONTRAST_TARGET = {}
    for r in contrasts.itertuples():
        for g in (r.group1, r.group2):
            if g in CONTROL_GROUPS:
                raise ValueError(
                    f"contrast {r.label}: '{g}' is a control group and cannot be contrasted."
                )
        t1, t2 = GROUP_INFO[r.group1]["target"], GROUP_INFO[r.group2]["target"]
        if t1 != t2:
            raise ValueError(
                f"contrast {r.label}: {r.group1} ({t1}) and {r.group2} ({t2}) are "
                "different targets; contrasts compare groups of one target."
            )
        CONTRAST_TARGET[r.label] = t1
    CONTRAST_GROUPS = {r.label: (r.group1, r.group2) for r in contrasts.itertuples()}
else:
    CONTRASTS, CONTRAST_TARGET, CONTRAST_GROUPS = [], {}, {}


def target_contrasts(target):
    return [c for c in CONTRASTS if CONTRAST_TARGET[c] == target]


# ---------------------------------------------------------------------------
# Wildcard constraints
# ---------------------------------------------------------------------------
def _alt(values):
    return "|".join(re.escape(v) for v in values) if values else "__none__"


CALLERS = ["macs2"] + (["seacr"] if PK["seacr"]["run"] else [])
if PK["caller"] not in CALLERS:
    raise ValueError(
        "peaks.caller is seacr, but peaks.seacr.run is false; switch SEACR on "
        "or use peaks.caller: macs2."
    )
PRIMARY = PK["caller"]
if PK["seacr"]["run"] and not (IGG["as_control"] and HAS_IGG):
    _warn(
        "SEACR without an IgG control uses a numeric threshold "
        f"(top {PK['seacr']['threshold']} of signal blocks); in the lab benchmark "
        "(mcf7 analysis 19) this gave degenerate peak sets. Prefer MACS2, or "
        "give SEACR an IgG control."
    )


wildcard_constraints:
    lib=_alt(LIBS),
    group=_alt(GROUPS),
    target=_alt(TARGETS),
    caller=_alt(CALLERS),
    level="individual|merged",
    name=_alt(LIBS + GROUPS),
    method="depth|greenlist|spikein",
    label=_alt(CONTRASTS),


# ---------------------------------------------------------------------------
# Reference helpers
# ---------------------------------------------------------------------------
REF = config["reference"]
BLACKLIST = REF.get("blacklist") or ""
USE_BLACKLIST = bool(BLACKLIST)
MITO = REF.get("mito_chrom", "chrM")
SPIKE = config.get("spikein") or {}
SPIKE_PREFIX = SPIKE.get("contig_prefix") or "spikein_"
SPIKE_ALIGN = INPUT_MODE == "fastq" and bool(
    SPIKE.get("fasta") or SPIKE.get("combined_bowtie2_index")
)
BUNDLED_GREENLISTS = {
    "hg38": "resources/greenlists/hg38_CUTnRUN_greenlist.v1.bed",
}


def greenlist_bed():
    g = REF.get("greenlist") or ""
    if g in BUNDLED_GREENLISTS:
        return os.path.join(workflow.basedir, BUNDLED_GREENLISTS[g])
    return g


def blacklist_input():
    return [BLACKLIST] if USE_BLACKLIST else []


def _index_files(prefix, detect):
    ext = ".bt2"
    if detect:
        found = glob.glob(f"{prefix}.1.bt2*")
        ext = ".bt2l" if found and found[0].endswith(".bt2l") else ".bt2"
    return [f"{prefix}.{s}{ext}" for s in ("1", "2", "3", "4", "rev.1", "rev.2")]


def bowtie2_index_files():
    """Index used for alignment: host, or host + spike-in when configured."""
    if SPIKE_ALIGN:
        if SPIKE.get("combined_bowtie2_index"):
            return _index_files(SPIKE["combined_bowtie2_index"], True)
        return _index_files("results/reference/bowtie2_combined/genome", False)
    if REF.get("bowtie2_index"):
        return _index_files(REF["bowtie2_index"], True)
    return _index_files("results/reference/bowtie2/genome", False)


# ---------------------------------------------------------------------------
# FASTQ-mode inputs
# ---------------------------------------------------------------------------
def _prefixed(p, base):
    if p and base and not os.path.isabs(p):
        return os.path.join(base, p)
    return p


def lib_runs(lib):
    rows = samples_raw[samples_raw["lib"] == lib]
    base = config.get("fastq_dir") or ""
    return [
        (_prefixed(r.fastq_1, base), _prefixed(r.fastq_2, base))
        for r in rows.itertuples()
    ]


def trim_inputs(wildcards):
    runs = lib_runs(wildcards.lib)
    if len(runs) == 1:
        return {"r1": runs[0][0], "r2": runs[0][1]}
    return {
        "r1": f"results/fastq/{wildcards.lib}_R1.merged.fastq.gz",
        "r2": f"results/fastq/{wildcards.lib}_R2.merged.fastq.gz",
    }


def merge_fastq_inputs(wildcards):
    idx = 0 if wildcards.read == "1" else 1
    return [r[idx] for r in lib_runs(wildcards.lib)]


def remove_dups(lib):
    D = config["duplicates"]
    return D["remove_control"] if is_control_lib(lib) else D["remove_target"]


# ---------------------------------------------------------------------------
# BAM-mode inputs
# ---------------------------------------------------------------------------
def import_bam_input(wildcards):
    return _prefixed(LIBRARIES.loc[wildcards.lib, "bam"], config.get("bam_dir") or "")


# ---------------------------------------------------------------------------
# BAM accessors (identical downstream of both input modes)
# ---------------------------------------------------------------------------
def lib_bam(lib):
    return f"results/bam/{lib}.final.bam"


def all_lib_bams(libs=None):
    return [lib_bam(l) for l in (libs if libs is not None else LIBS)]


PEAK_MAX_SIZE = int(config["fragments"].get("peak_max_size") or 0)


def size_selected(group):
    """TF groups are called on fragments <= fragments.peak_max_size when set."""
    return PEAK_MAX_SIZE > 0 and GROUP_INFO[group]["target_type"] == "tf"


def peak_bam(name):
    """BAM a peak set is called on: a library, or a group's merged BAM."""
    if name in LIBRARIES.index:
        if size_selected(lib_group(name)):
            return f"results/bam/sized/{name}.bam"
        return lib_bam(name)
    kind = "merged_sized" if size_selected(name) else "merged"
    return f"results/bam/{kind}/{name}.bam"


def name_group(name):
    return lib_group(name) if name in LIBRARIES.index else name


def _merged_kind(name):
    return "merged_sized" if size_selected(name_group(name)) else "merged"


def control_bam(name):
    """Merged control BAM for MACS2 -c / SEACR, or [] (no control)."""
    c = group_control(name_group(name))
    return [f"results/bam/{_merged_kind(name)}/{c}.bam"] if c else []


def igg_reference_bam(name):
    """IgG BAM for the enrichment gate: the group's control, else all IgG pooled."""
    c = GROUP_INFO[name_group(name)]["control"]
    if c:
        return f"results/bam/{_merged_kind(name)}/{c}.bam"
    sfx = "_sized" if size_selected(name_group(name)) else ""
    return f"results/bam/igg_pool/pool{sfx}.bam"


def merge_group_inputs(wildcards):
    libs = LIBS_BY_GROUP[wildcards.group]
    if wildcards.kind == "merged_sized":
        return [f"results/bam/sized/{l}.bam" for l in libs]
    return all_lib_bams(libs)


def igg_pool_inputs(wildcards):
    if wildcards.sfx == "_sized":
        return [f"results/bam/sized/{l}.bam" for l in CONTROL_LIBS]
    return all_lib_bams(CONTROL_LIBS)


# ---------------------------------------------------------------------------
# Peaks
# ---------------------------------------------------------------------------
def peak_names(level):
    return TARGET_LIBS if level == "individual" else TARGET_GROUPS


def macs2_mode_args(name, relaxed=False):
    g = name_group(name)
    if GROUP_INFO[g]["peak_mode"] == "broad":
        return f"--broad --broad-cutoff {PK['broad_cutoff']}"
    if relaxed:
        return f"-p {PK['relaxed_pvalue']}"
    return f"-q {PK['qvalue']}"


def peaks_final(caller, level, name):
    return f"results/peaks/{caller}/{level}/{name}.bed"


def group_rep_threshold(group):
    return min(int(PK["min_overlap"]), len(LIBS_BY_GROUP[group]))


SINGLE_REP_GROUPS = [g for g in TARGET_GROUPS if len(LIBS_BY_GROUP[g]) < 2]
if SINGLE_REP_GROUPS:
    _warn(
        "single-replicate group(s): "
        + ", ".join(SINGLE_REP_GROUPS)
        + "; their consensus peaks are the one library's peaks (no reproducibility)."
    )

# ---------------------------------------------------------------------------
# IDR (MACS2, narrow groups with >= 2 replicates)
# ---------------------------------------------------------------------------
IDR = config["idr"]
IDR_GROUPS = [
    g
    for g in TARGET_GROUPS
    if len(LIBS_BY_GROUP[g]) >= 2 and GROUP_INFO[g]["peak_mode"] == "narrow"
]


def idr_pairs(group):
    libs = LIBS_BY_GROUP[group]
    if len(libs) > 2 and IDR.get("multi_rep_strategy") == "first_two":
        return [(libs[0], libs[1])]
    return list(itertools.combinations(libs, 2))


def idr_pair_files(wildcards):
    return [
        f"results/peaks/idr/pairs/{wildcards.group}/{a}__{b}.idr.narrowPeak"
        for a, b in idr_pairs(wildcards.group)
    ]


def idr_scaled_threshold():
    # IDR column 5 = min(int(-125 * log2(IDR)), 1000); 0.05 -> 540
    import math

    return min(int(-125 * math.log2(float(IDR["threshold"]))), 1000)


def idr_group_peaks(group):
    """IDR set for IDR groups; merged peaks for single-rep or broad groups."""
    if group in IDR_GROUPS:
        return f"results/peaks/idr/{group}_idr.narrowPeak"
    return peaks_final("macs2", "merged", group)


# ---------------------------------------------------------------------------
# Peak sets used by the analyses (per target)
# ---------------------------------------------------------------------------
PEAKSET_KINDS = {
    "consensus": "reproducible peaks (>= peaks.min_overlap replicates) per group, merged over the target's groups",
    "merged": "peaks called on each group's merged BAM, merged over the target's groups",
    "idr": "IDR-reproducible peaks per group (modules.idr), merged over the target's groups",
}
DIFFBIND_INDIVIDUAL = "individual"


def target_peakset(kind, target):
    return f"results/peaks/consensus/{target}.{kind}.bed"


def _analysis_on(section):
    return bool(config[section].get("run", False))


def resolve_peaks(section, allow_individual=False):
    """(label, target -> BED or None) for `<section>.peaks`."""
    value = str(config[section].get("peaks") or "").strip()
    choices = list(PEAKSET_KINDS) + ([DIFFBIND_INDIVIDUAL] if allow_individual else [])
    hint = f"one of {', '.join(choices)}, or a path to a BED file"
    if not value:
        raise ValueError(
            f"{section}.run is true, so {section}.peaks must be set ({hint})."
        )
    if value == DIFFBIND_INDIVIDUAL:
        if not allow_individual:
            raise ValueError(
                f"{section}.peaks: '{DIFFBIND_INDIVIDUAL}' is only valid for diff ({hint})."
            )
        return value, lambda t: None
    if value in PEAKSET_KINDS:
        if value == "idr" and not MODULES.get("idr", False):
            raise ValueError(f"{section}.peaks is idr, but modules.idr is false.")
        return value, lambda t: target_peakset(value, t)
    if not os.path.isfile(value):
        raise ValueError(
            f"{section}.peaks: '{value}' is neither a pipeline peak set nor an "
            f"existing file ({hint})."
        )
    stem = re.sub(r"[^A-Za-z0-9]+", "_", Path(value).name.split(".")[0]).strip("_")
    return f"custom_{stem}", lambda t: value


RUN_DIFF = _analysis_on("diff")
RUN_NORMCHECK = _analysis_on("normcheck")
RUN_HEATMAPS = _analysis_on("heatmaps")
RUN_ANNOTATE = _analysis_on("annotate")
RUN_MOTIFS = _analysis_on("motifs")

if (RUN_DIFF or RUN_NORMCHECK) and not CONTRASTS:
    raise ValueError(
        "diff/normcheck are on, but no contrasts are defined (config key `contrasts`)."
    )
if RUN_NORMCHECK and not RUN_DIFF:
    raise ValueError(
        "normcheck.run needs diff.run: it re-runs the diff contrasts on the same "
        "counted DiffBind object under other normalizations."
    )
if RUN_DIFF:
    for _lab, (_g1, _g2) in CONTRAST_GROUPS.items():
        for _g in (_g1, _g2):
            if len(LIBS_BY_GROUP[_g]) < 2:
                raise ValueError(
                    f"contrast {_lab}: group '{_g}' has one replicate; DiffBind "
                    "(DESeq2) needs >= 2 replicates per contrasted group."
                )
if RUN_ANNOTATE and not REF.get("gtf"):
    raise ValueError("annotate.run needs reference.gtf (the TxDb is built from it).")

DIFF_PEAKS = resolve_peaks("diff", allow_individual=True) if RUN_DIFF else None
HM_PEAKS = resolve_peaks("heatmaps") if RUN_HEATMAPS else None
ANN_PEAKS = resolve_peaks("annotate") if RUN_ANNOTATE else None
MOTIF_PEAKS = resolve_peaks("motifs") if RUN_MOTIFS else None
DIFF_TARGETS = (
    list(dict.fromkeys(CONTRAST_TARGET[c] for c in CONTRASTS)) if RUN_DIFF else []
)

# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------
DIFF = config["diff"]
NC = config["normcheck"]
DIFF_NORM = DIFF.get("normalization", "depth")
NC_METHODS = (
    [m for m in NC.get("methods", []) if m != DIFF_NORM] if RUN_NORMCHECK else []
)

NORM_METHODS_NEEDED = set(NORM["bigwig_methods"])
if RUN_DIFF:
    NORM_METHODS_NEEDED.add(DIFF_NORM)
NORM_METHODS_NEEDED |= {m for m in NC_METHODS if m != "csaw"}
if RUN_HEATMAPS:
    NORM_METHODS_NEEDED.add(config["heatmaps"].get("normalization", "depth"))
SIZEFACTOR_METHODS = sorted(NORM_METHODS_NEEDED - {"depth"})

HAS_SPIKEIN = SPIKE_ALIGN or (
    INPUT_MODE == "bam" and bool((LIBRARIES["spikein_reads"] != "").all())
)
if "spikein" in NORM_METHODS_NEEDED and not HAS_SPIKEIN:
    raise ValueError(
        "spike-in normalization is requested (normalization.bigwig_methods, "
        "diff.normalization, normcheck.methods or heatmaps.normalization), but "
        "there are no spike-in counts: set spikein.fasta or "
        "spikein.combined_bowtie2_index (fastq mode), or fill the samplesheet "
        "`spikein_reads` column for every library (bam mode)."
    )
if "greenlist" in NORM_METHODS_NEEDED and not greenlist_bed():
    raise ValueError(
        "greenlist normalization is requested, but reference.greenlist is empty. "
        f"Use a bundled greenlist ({', '.join(BUNDLED_GREENLISTS)}) or a BED path."
    )


def norm_group(lib):
    t = LIBRARIES.loc[lib, "target"]
    if NORM.get("greenlist_group") == "target_batch":
        return f"{t}.{LIBRARIES.loc[lib, 'batch'] or 'nobatch'}"
    return t


def bigwig(method, lib):
    return f"results/bigwig/{method}/{lib}.bw"


def group_bigwig(method, group):
    return f"results/bigwig/{method}/groups/{group}.bw"


def sizefactor_table(method):
    return f"results/normalization/{method}/size_factors.tsv"


def spikein_count_files():
    if SPIKE_ALIGN:
        return [f"results/qc/spikein/{l}.spikein.tsv" for l in LIBS]
    return []


# ---------------------------------------------------------------------------
# Batch covariate (config `batch: true` + samplesheet `batch` column)
# ---------------------------------------------------------------------------
USE_BATCH = bool(config.get("batch", False))


def _check_design(label, libs):
    """Stop unless ~batch + group over `libs` is full rank with residual df."""
    import numpy as np

    meta = LIBRARIES.loc[libs]
    x = pd.get_dummies(meta[["group", "batch"]], drop_first=True, dtype=float)
    x.insert(0, "intercept", 1.0)
    rank = np.linalg.matrix_rank(x.to_numpy())
    if rank < x.shape[1]:
        raise ValueError(
            f"batch: {label}: ~batch + group is not estimable, batch is "
            "confounded with group. Libraries per group and batch:\n"
            + pd.crosstab(meta["group"], meta["batch"]).to_string()
        )
    if len(libs) <= rank:
        raise ValueError(
            f"batch: {label}: {len(libs)} libraries leave no residual degrees "
            f"of freedom for ~batch + group ({rank} coefficients)."
        )


def diff_libs(target):
    return [l for g in GROUPS_BY_TARGET[target] for l in LIBS_BY_GROUP[g]]


if USE_BATCH and RUN_DIFF:
    for _t in DIFF_TARGETS:
        _libs = diff_libs(_t)
        _nobatch = [l for l in _libs if not LIBRARIES.loc[l, "batch"]]
        if _nobatch:
            raise ValueError(
                "batch: true needs a `batch` value for every library in the "
                f"statistics; missing for: {', '.join(_nobatch)}"
            )
        if LIBRARIES.loc[_libs, "batch"].nunique() > 1:
            _check_design(f"target {_t}", _libs)

# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------
# CUT&RUN thresholds: lab STANDARDS.md (aligned > 85 %, mito < 5 %,
# duplicates < 50 %, TF libraries mostly < 120 bp fragments); fragments and
# FRiP are CUT&RUN-typical (Meers et al. 2019; nf-core/cutandrun).
QC_THRESHOLDS = {
    "fragments": [5000000, 2000000],
    "aligned_pct": [85, 70],
    "mito_pct": [5, 10],
    "dup_pct": [50, 70],
    "frip": [0.1, 0.05],
    "tf_fraction": [0.5, 0.3],
    "igg_pass_fraction": [0.8, 0.5],
}
QC_THRESHOLDS.update(config["qc"].get("thresholds") or {})


def frip_peaks(lib):
    if is_control_lib(lib):
        return []
    return peaks_final(PRIMARY, "merged", lib_group(lib))


def qc_summary_inputs(wildcards):
    files = {
        "flagstat": expand("results/qc/flagstat/{lib}.flagstat.txt", lib=LIBS),
        "frip": expand("results/qc/frip/{lib}.frip.tsv", lib=TARGET_LIBS),
        "frag": expand("results/qc/fragment_sizes/{lib}_fragment_sizes.tsv", lib=LIBS),
    }
    if INPUT_MODE == "fastq":
        files["report"] = "results/qc/alignment_qc_report.tsv"
    if IGG_FOLD:
        files["igg"] = [
            f"results/peaks/{PRIMARY}/individual/{l}.igg_summary.tsv"
            for l in TARGET_LIBS
        ]
    if HAS_SPIKEIN:
        files["spikein"] = "results/qc/spikein_summary.tsv"
    return files


def multiqc_inputs(wildcards):
    files = expand("results/qc/flagstat/{lib}.flagstat.txt", lib=LIBS)
    if INPUT_MODE == "fastq":
        files += expand("logs/align/{lib}.cutadapt.log", lib=LIBS)
        files += expand("logs/align/{lib}.bowtie2.log", lib=LIBS)
        files += expand("results/qc/markdup/{lib}.markdup.txt", lib=LIBS)
        files.append("results/qc/alignment_qc_report.tsv")
    files.append("results/qc/qc_summary.tsv")
    if MODULES.get("qc", True):
        if INPUT_MODE == "fastq":
            files += expand(
                "results/qc/fastqc/{lib}_R{r}_fastqc.zip", lib=LIBS, r=["1", "2"]
            )
        files += [
            "results/qc/deeptools/fragment_size_table.tsv",
            "results/qc/deeptools/fragment_size_raw.tsv",
            "results/qc/deeptools/correlation_spearman.tsv",
            "results/qc/deeptools/pca_data.tsv",
            "results/qc/deeptools/fingerprint_metrics.tsv",
            "results/qc/deeptools/fingerprint_counts.tsv",
        ]
    files += expand("results/peaks/macs2/individual/{lib}_peaks.xls", lib=TARGET_LIBS)
    for caller in CALLERS:
        files += [
            f"results/peaks/{caller}/individual/peak_summary.tsv",
            f"results/peaks/{caller}/merged/peak_summary.tsv",
        ]
    if IGG_FOLD:
        files.append("results/qc/igg_enrichment.tsv")
    if HAS_SPIKEIN:
        files.append("results/qc/spikein_summary.tsv")
    for m in SIZEFACTOR_METHODS:
        files.append(sizefactor_table(m))
    if MODULES.get("idr", False) and IDR_GROUPS:
        files.append("results/peaks/idr/idr_summary.tsv")
    return files


# ---------------------------------------------------------------------------
# Rule input helpers (peaks, downstream)
# ---------------------------------------------------------------------------
def target_set_inputs(wildcards):
    groups = GROUPS_BY_TARGET[wildcards.target]
    if wildcards.kind == "consensus":
        return [f"results/peaks/consensus/groups/{g}.bed" for g in groups]
    if wildcards.kind == "merged":
        return [peaks_final(PRIMARY, "merged", g) for g in groups]
    return [idr_group_peaks(g) for g in groups]


def heatmap_groups(target):
    groups = list(GROUPS_BY_TARGET[target])
    controls = [GROUP_INFO[g]["control"] for g in groups if GROUP_INFO[g]["control"]]
    return groups + list(dict.fromkeys(controls))


def heatmap_regions(wildcards):
    if wildcards.kind == "peaks":
        return [HM_PEAKS[1](wildcards.target)] if HM_PEAKS else []
    return ["results/reference/tss.bed"]


def contrast_sets():
    if not RUN_DIFF:
        return {}
    return {
        f"{c}_{d}": f"results/diff/contrasts/{c}_{d}.bed"
        for c in CONTRASTS
        for d in ("gained", "lost")
    }


# ---------------------------------------------------------------------------
# Resources
# ---------------------------------------------------------------------------
def threads(key, default=4):
    return int(config.get("threads", {}).get(key, default))


# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------
def core_targets():
    t = all_lib_bams() + [f"{b}.bai" for b in all_lib_bams()]
    if INPUT_MODE == "fastq":
        t.append("results/qc/alignment_qc_report.tsv")
    t += [f"results/qc/fragment_sizes/{l}_fragment_sizes.tsv" for l in LIBS]
    for caller in CALLERS:
        t += [
            f"results/peaks/{caller}/individual/peak_summary.tsv",
            f"results/peaks/{caller}/merged/peak_summary.tsv",
        ]
    t += [target_peakset(k, tg) for tg in TARGETS for k in ("consensus", "merged")]
    t.append("results/peaks/consensus/consensus_summary.tsv")
    if MODULES.get("idr", False):
        t += [target_peakset("idr", tg) for tg in TARGETS]
        if IDR_GROUPS:
            t.append("results/peaks/idr/idr_summary.tsv")
    if IGG_FOLD:
        t.append("results/qc/igg_enrichment.tsv")
    if HAS_SPIKEIN:
        t.append("results/qc/spikein_summary.tsv")
    return t


def processing_targets():
    """Default target: processed data (BAMs, bigWigs, peaks, QC)."""
    t = core_targets()
    t.append("results/qc/qc_summary.tsv")
    for m in SIZEFACTOR_METHODS:
        t.append(sizefactor_table(m))
    for m in NORM["bigwig_methods"]:
        t += [bigwig(m, l) for l in LIBS]
        t += [group_bigwig(m, g) for g in GROUPS]
    if MODULES.get("qc", True):
        t += [
            "results/qc/deeptools/fragment_size_distribution.png",
            "results/qc/deeptools/correlation_spearman.png",
            "results/qc/deeptools/pca_plot.png",
            "results/qc/deeptools/fingerprint.png",
        ]
    if MODULES.get("multiqc", True):
        t.append("results/qc/multiqc/multiqc_report.html")
    return t


def analysis_targets():
    """Outputs of the analyses switched on with `<section>.run: true`."""
    t = []
    for tg in DIFF_TARGETS:
        t.append(f"results/diff/{tg}/{DIFF_NORM}/summary.tsv")
        if RUN_NORMCHECK:
            t.append(f"results/diff/{tg}/normcheck/norm_verdict.tsv")
    if RUN_HEATMAPS:
        t += [f"results/heatmaps/{tg}/peaks_heatmap.png" for tg in TARGETS]
        if config["heatmaps"].get("tss", True) and REF.get("gtf"):
            t += [f"results/heatmaps/{tg}/tss_heatmap.png" for tg in TARGETS]
        if RUN_DIFF:
            t += [f"results/heatmaps/contrasts/{c}_heatmap.png" for c in CONTRASTS]
    if RUN_ANNOTATE:
        t.append("results/annotate/annotation_summary.tsv")
    if RUN_MOTIFS:
        t += [f"results/motifs/{tg}/knownResults.txt" for tg in TARGETS]
        if RUN_DIFF:
            t += [
                f"results/motifs/contrasts/{c}_{d}/knownResults.txt"
                for c in CONTRASTS
                for d in ("gained", "lost")
            ]
    return t


def all_targets():
    return processing_targets() + analysis_targets()
