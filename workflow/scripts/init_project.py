#!/usr/bin/env python3
"""Set up a CUT&RUN project config from the lab's reference manifest.

Asks which genome, whether there are IgG controls, which spike-in and which
blacklist, then writes <project>/project.yaml with the matching paths from the
resource manifest (/data/resource/manifest.yaml on the Salk server), ready for

    pixi run snakemake -s workflow/Snakefile --directory <project> \
        --configfile <project>/project.yaml --profile profiles/local -n

Run it interactively (`pixi run init`), or give every answer as a flag, which is
how agents should call it after asking the user:

    pixi run init --dir /path/to/project --genome hg38 --igg yes \
        --spikein ecoli --blacklist encode_v2

`--list` prints the choices the manifest offers. The manifest is found at
--manifest, else $MCBLA_RESOURCE_MANIFEST, else /data/resource/manifest.yaml.
"""

import argparse
import datetime
import json
import os
import shutil
import sys

import yaml

DEFAULT_MANIFEST = "/data/resource/manifest.yaml"
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
# pixi runs tasks from the repo root; INIT_CWD is where the user typed the command
CWD = os.environ.get("INIT_CWD") or os.getcwd()


def die(msg, code=1):
    print(f"init: {msg}", file=sys.stderr)
    sys.exit(code)


def find_manifest(arg):
    for p in (arg, os.environ.get("MCBLA_RESOURCE_MANIFEST"), DEFAULT_MANIFEST):
        if p and os.path.exists(p):
            return p
    die(
        "no resource manifest found (looked at --manifest, $MCBLA_RESOURCE_MANIFEST, "
        f"{DEFAULT_MANIFEST}). Outside the Salk server, copy config/salk_example.yaml "
        "and fill in your reference paths by hand (see README)."
    )


def usable(entry, allow_unverified):
    return entry.get("status", "ok") == "ok" or allow_unverified


def choices(genome, allow_unverified):
    """Blacklist and spike-in keys offered for one genome."""
    bl = {
        k: v
        for k, v in (genome.get("blacklists") or {}).items()
        if usable(v, allow_unverified)
    }
    sp = {
        k: v
        for k, v in (genome.get("spikeins") or {}).items()
        if usable(v, allow_unverified)
    }
    return bl, sp


def default_blacklist(bl):
    for k, v in bl.items():
        if v.get("default"):
            return k
    return next(iter(bl), "none")


def resolve_genome(genomes, answer):
    a = answer.lower()
    for key, g in genomes.items():
        if a == key.lower() or a in (x.lower() for x in g.get("aliases", [])):
            return key
    die(f"unknown genome '{answer}'; choose one of: {', '.join(genomes)}")


def ask(question, options, default=None):
    """Numbered-choice prompt. options: list of (value, description)."""
    print(f"\n{question}")
    for i, (val, desc) in enumerate(options, 1):
        mark = "  (default)" if val == default else ""
        print(f"  {i}. {val:<16} {desc}{mark}")
    while True:
        raw = input(
            f"Choice [1-{len(options)}{', Enter = default' if default else ''}]: "
        ).strip()
        if not raw and default:
            return default
        if raw.isdigit() and 1 <= int(raw) <= len(options):
            return options[int(raw) - 1][0]
        for val, _ in options:
            if raw.lower() == val.lower():
                return val
        if os.path.sep in raw or raw.endswith(".bed"):  # custom BED path
            return raw
        print("  Please enter one of the numbers above.")


def ask_text(question, default):
    raw = input(f"\n{question} [{default}]: ").strip()
    return raw or default


def list_choices(man, genomes, allow_unverified, as_json):
    out = {}
    for key, g in genomes.items():
        bl, sp = choices(g, allow_unverified)
        out[key] = {
            "description": g.get("description", ""),
            "blacklists": {k: v.get("description", "") for k, v in bl.items()},
            "default_blacklist": default_blacklist(bl),
            "spikeins": {k: v.get("description", "") for k, v in sp.items()},
            "greenlist": bool(g.get("greenlist")),
        }
    if as_json:
        print(json.dumps(out, indent=2))
        return
    print(
        f"Manifest version {man['manifest_version']} (updated {man.get('updated')})\n"
    )
    for key, o in out.items():
        print(f"{key}: {o['description']}")
        bls = [
            f"{k} (default)" if k == o["default_blacklist"] else k
            for k in o["blacklists"]
        ]
        print("  blacklists: " + ", ".join(bls + ["none"]))
        print("  spike-ins:  " + ", ".join(list(o["spikeins"]) + ["none"]))
        print(f"  greenlist:  {'yes' if o['greenlist'] else 'no'}\n")


def check_exists(paths):
    missing = [p for p in paths if p and not os.path.exists(p)]
    if missing:
        die(
            "these resources are missing on disk (tell the manifest maintainer):\n  "
            + "\n  ".join(missing)
        )


def bt2(prefix):
    for ext in (".bt2", ".bt2l"):
        if os.path.exists(prefix + ".1" + ext):
            return prefix + ".1" + ext
    return prefix + ".1.bt2"


def q(v):
    """YAML scalar."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    return json.dumps(str(v))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--dir", help="project directory (created if needed)")
    ap.add_argument(
        "--genome", help="genome key or alias from the manifest (hg38, mm10, mm39, ...)"
    )
    ap.add_argument(
        "--igg", choices=["yes", "no"], help="are there IgG control libraries?"
    )
    ap.add_argument(
        "--spikein", help="spike-in key from the manifest (ecoli, ...) or none"
    )
    ap.add_argument(
        "--blacklist", help="blacklist key from the manifest, none, or a BED path"
    )
    ap.add_argument(
        "--fastq-dir", help="directory holding the FASTQs (default: <project>/fastq)"
    )
    ap.add_argument(
        "--manifest", help=f"resource manifest (default {DEFAULT_MANIFEST})"
    )
    ap.add_argument(
        "--allow-unverified",
        action="store_true",
        help="also offer entries not marked status: ok",
    )
    ap.add_argument(
        "--force", action="store_true", help="overwrite an existing project.yaml"
    )
    ap.add_argument(
        "--list", action="store_true", help="print the manifest's choices and exit"
    )
    ap.add_argument(
        "--json", action="store_true", help="with --list: machine-readable output"
    )
    a = ap.parse_args()

    mpath = find_manifest(a.manifest)
    with open(mpath) as fh:
        man = yaml.safe_load(fh)
    genomes = {
        k: g
        for k, g in man["genomes"].items()
        if g.get("bowtie2_index") and usable(g, a.allow_unverified)
    }
    if a.list:
        list_choices(man, genomes, a.allow_unverified, a.json)
        return

    interactive = sys.stdin.isatty()
    needed = [
        f
        for f, v in (
            ("--dir", a.dir),
            ("--genome", a.genome),
            ("--igg", a.igg),
            ("--spikein", a.spikein),
            ("--blacklist", a.blacklist),
        )
        if v is None
    ]
    if needed and not interactive:
        die(
            "not a terminal, so every answer must be a flag. Missing: "
            + ", ".join(needed)
            + ".\n"
            "Ask the user (genome; IgG controls yes/no; spike-in; blacklist) -- do not guess. "
            "`pixi run init --list` shows the options.",
            2,
        )

    if a.dir is None:
        a.dir = ask_text("Project directory (created if needed)", ".")
    proj = os.path.abspath(os.path.join(CWD, a.dir))
    out = os.path.join(proj, "project.yaml")
    if os.path.exists(out) and not a.force:
        die(f"{out} exists; use --force to overwrite")

    if a.genome is None:
        a.genome = ask(
            "Which genome are the samples from?",
            [(k, g.get("description", "")) for k, g in genomes.items()],
            "hg38" if "hg38" in genomes else None,
        )
    gkey = resolve_genome(genomes, a.genome)
    g = genomes[gkey]
    bl, sp = choices(g, a.allow_unverified)

    if a.igg is None:
        a.igg = ask(
            "Do you have IgG control libraries?",
            [
                ("yes", "call peaks against IgG and report enrichment over IgG"),
                ("no", "call peaks without a control"),
            ],
            "yes",
        )

    if a.spikein is None:
        a.spikein = ask(
            "Spike-in for normalization?",
            [(k, v.get("description", "")) for k, v in sp.items()]
            + [("none", "no spike-in (depth / greenlist normalization only)")],
            "none",
        )
    if a.spikein != "none" and a.spikein not in sp:
        die(
            f"spike-in '{a.spikein}' is not available for {gkey}; choose: {', '.join(list(sp) + ['none'])}"
        )

    if a.blacklist is None:
        a.blacklist = ask(
            "Blacklist to filter reads and peaks? (or type a BED path)",
            [(k, v.get("description", "")) for k, v in bl.items()]
            + [("none", "no blacklist filtering")],
            default_blacklist(bl),
        )
    if a.blacklist in bl:
        blacklist = bl[a.blacklist]["path"]
    elif a.blacklist == "none":
        blacklist = ""
    elif os.path.exists(os.path.join(CWD, a.blacklist)):
        blacklist = os.path.abspath(os.path.join(CWD, a.blacklist))
    else:
        die(
            f"blacklist '{a.blacklist}' is neither a {gkey} manifest key ({', '.join(bl) or 'none'}), 'none', nor an existing file"
        )

    if a.fastq_dir is None:
        a.fastq_dir = (
            ask_text(
                "Directory with the FASTQs (relative to the project directory, or absolute)",
                "fastq",
            )
            if interactive
            else "fastq"
        )

    greenlist = (g.get("greenlist") or {}).get("path", "")
    spike = sp.get(a.spikein) if a.spikein != "none" else None
    check_exists(
        [
            g["fasta"],
            g["fasta"] + ".fai",
            bt2(g["bowtie2_index"]),
            g.get("gtf"),
            blacklist,
            greenlist,
            bt2(spike["combined_bowtie2_index"]) if spike else None,
        ]
    )

    methods = ["depth"] + (["spikein"] if spike else [])
    igg = a.igg == "yes"
    run = (
        f"pixi run snakemake -s workflow/Snakefile --directory {proj} "
        f"--configfile {out} --profile profiles/local -n"
    )
    lines = [
        "# =============================================================================",
        f"# CUT&RUN project config -- written by `pixi run init` on {datetime.date.today()}",
        f"# from {mpath} (manifest version {man['manifest_version']}, updated {man.get('updated')})",
        f"# Choices: genome={gkey}  igg={a.igg}  spikein={a.spikein}  blacklist={a.blacklist}",
        "#",
        "# Holds only the keys that differ from the pipeline's config/config.yaml.",
        "# Relative paths are relative to this directory. Run from the pipeline repo:",
        f"#   {run}",
        "# =============================================================================",
        "",
        "# Fill these in (see the pipeline README, 'Samplesheet').",
        "samples: samples.tsv",
        "contrasts: contrasts.tsv",
        f"fastq_dir: {q(a.fastq_dir)}",
        "",
        f"# {g.get('description', gkey)}",
        "reference:",
        f"  fasta: {q(g['fasta'])}",
        f"  bowtie2_index: {q(g['bowtie2_index'])}",
        f"  gtf: {q(g.get('gtf') or '')}",
        f"  blacklist: {q(blacklist)}"
        + (f"    # {a.blacklist}" if a.blacklist in bl else ""),
        f"  greenlist: {q(greenlist)}"
        + ("" if greenlist else "    # no CUT&RUN greenlist for this genome"),
        f"  effective_genome_size: {q(g['effective_genome_size'])}",
        f"  macs2_gsize: {q(g['macs2_gsize'])}",
        f"  mito_chrom: {q(g.get('mito_chrom', 'chrM'))}",
        "",
    ]
    if spike:
        lines += [
            f"# {spike.get('description', a.spikein)}",
            "spikein:",
            f"  combined_bowtie2_index: {q(spike['combined_bowtie2_index'])}",
            f"  contig_prefix: {q(spike['contig_prefix'])}",
            "",
        ]
    lines += [
        "# IgG: "
        + (
            "fill the samplesheet `control` column with the IgG group's name"
            if igg
            else "no IgG libraries; peaks are called without a control"
        ),
        "igg:",
        f"  as_control: {q(igg)}",
        f"  qc: {q(igg)}",
        "",
        "normalization:",
        f"  bigwig_methods: [{', '.join(methods)}]"
        + (
            ""
            if not greenlist
            else "    # add greenlist for contrasts with a global shift"
        ),
        "",
    ]

    os.makedirs(proj, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("\n".join(lines))
    copied = []
    for name in ("samples.tsv", "contrasts.tsv"):
        dst = os.path.join(proj, name)
        if not os.path.exists(dst):
            shutil.copy(os.path.join(REPO, "config", name), dst)
            copied.append(dst)

    print(f"\nWrote {out}")
    for c in copied:
        what = "libraries" if c.endswith("samples.tsv") else "contrasts"
        print(
            f"Copied the example {os.path.basename(c)} to {c} -- replace its rows with your {what}."
        )
    if igg:
        print(
            "IgG: set the `control` column of every target row to the IgG group's name."
        )
    if spike:
        print(
            "Spike-in: add `spikein` to diff.normalization in project.yaml to use it for DiffBind."
        )
    print(f"\nNext, a dry run from {REPO}:\n  {run}")


if __name__ == "__main__":
    main()
