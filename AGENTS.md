# AGENTS.md — running mcbla-cutrun-pipe for a lab member

Instructions for AI agents (Claude Code, Codex, …) asked to run this CUT&RUN
pipeline on the McBla lab server. People can follow them too. For pipeline
internals, read `README.md` and `docs/`.

## 1. Before running anything, ask the user

Ask these questions and wait for the answers. **Do not guess or pick defaults
without asking**, even if the FASTQ names seem to hint at an answer.

1. **Which genome?** Run `pixi run init --list` and show the options
   (`hg38`, `mm10`, `mm39`, `rn6`). Mention that mm10 and mm39 are different
   builds, so the user should confirm which one their other data uses.
2. **Are there IgG control libraries?** If yes, which samplesheet group holds
   them, and which target groups does each IgG control?
3. **Spike-in?** Choose from the options `--list` gives for that genome
   (`ecoli` = E. coli carry-over DNA from pA/pAG-MNase), or `none`.
4. **Blacklist?** The genome's default (ENCODE v2 for hg38/mm10, excluderanges
   for mm39), a CUT&RUN-specific list where one exists (`cutrun_demello`),
   `none`, or a path to their own BED. rn6 has no blacklist.

Also ask where the FASTQs are, where the project directory should go (usually
under their own `/data/<user>/`, never inside this repo or `/data/resource`),
and which groups to compare if they want differential binding.

## 2. Write the project config with `init`

Pass the user's answers as flags:

```bash
cd /path/to/mcbla-cutrun-pipe
pixi run init --dir /data/<user>/<project> --genome hg38 --igg yes \
    --spikein ecoli --blacklist encode_v2 --fastq-dir /path/to/fastqs
```

`init` takes the reference paths from the lab manifest,
`/data/resource/manifest.yaml`. It checks that each file exists, then writes
`<project>/project.yaml` and copies example `samples.tsv` / `contrasts.tsv`
into the project. It exits with status 2 if an answer is missing; that means
go back and ask the user.

- **Never download or build a genome, index or blacklist yourself**, and never
  point the config at files in someone's personal directory. If something the
  user needs isn't in the manifest (another genome, a dm6 or yeast spike-in),
  stop and tell them. New resources go into `/data/resource` through its
  maintainer (see `/data/resource/_admin/docs/README.md`).
- Entries marked `unverified` in the manifest are hidden. Only use
  `--allow-unverified` if the user explicitly asks for it after you explain
  why the entry is unverified.

## 3. Fill in the samplesheet

Edit `<project>/samples.tsv` (columns are described in the README,
"Samplesheet"). Show the finished sheet to the user before running. In
particular:

- With IgG: set `control` on every target row to the IgG group's name, and set
  `target` to `IgG` on the IgG rows.
- Set `target_type` to `tf` or `histone`. Set `peak_mode` to `broad` for broad
  marks (H3K27me3, H3K36me3, H3K9me3).
- Fill in `contrasts.tsv` only if `diff`/`normcheck` are switched on. Both
  groups in a contrast must have the same target and at least 2 replicates.

## 4. Dry run, then run

```bash
pixi run snakemake -s workflow/Snakefile --directory /data/<user>/<project> \
    --configfile /data/<user>/<project>/project.yaml -n            # dry run: fix any errors first
pixi run snakemake -s workflow/Snakefile --directory /data/<user>/<project> \
    --configfile /data/<user>/<project>/project.yaml --profile profiles/local --cores 32
```

- This server has **no Slurm**, so use `profiles/local`. It is shared, so keep
  `--cores` at 32 or fewer unless the user says otherwise.
- The run takes hours. Start it in the background or under `tmux`/`nohup`, and
  tell the user where the log is (`<project>/.snakemake/log/`).
- If the dry run warns about IgG from another batch or SEACR without IgG, pass
  the warning on to the user instead of silently changing settings.

## 5. Things that need the user's decision

- **Normalization:** `init` writes depth tracks, plus spike-in tracks if a
  spike-in was chosen. Greenlist normalization (hg38, mm39) and
  `diff.normalization` are choices for the user. Explain them using the
  README's "Normalization" section and don't switch them on unasked.
- **Downstream analyses** (`diff`, `normcheck`, `heatmaps`, `annotate`,
  `motifs`) are off by default. Turn one on only when the user asks for it.
