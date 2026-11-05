# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Snakemake pipeline for paired-end bulk RNA-seq (Daniel Fischer, Natural Resources Institute
Finland): per-lane FastQC, lane concatenation, adapter/quality trimming (fastp), optional
sequential decontamination against configurable reference genomes (STAR), genome+transcriptome
alignment (STAR, with a metatranscriptome mode that skips the annotation), transcript assembly
and cross-sample merging (StringTie), gene/transcript quantification (featureCounts against
both the reference annotation and the merged StringTie assembly; RSEM available on demand), and
one consolidated HTML report. It follows the same module/escalation/docker/report architecture
as `Snakebite-Holoruminant-MetaG` and `Snakebite-miRNA` (sibling pipelines in the same
"Snakebite" family) - those repos' `CLAUDE.md` files are useful cross-references for the shared
conventions. It runs on HPC (SLURM) with every tool delivered as a Singularity/Apptainer
container.

## Running the pipeline

There is no build step, and the only test is `tests/resources/check_resource_sets.py`, which
checks the resource sets against the published limits of CSC's Roihu. All "development" happens by
editing `.smk` rule files and dry-running Snakemake.

The pipeline is invoked through the wrapper script, never by calling `snakemake` directly:

```bash
bash run_Snakebite-RNAseq.sh <module_or_rule> [snakemake args...]

# dry run (always do this first when editing rules)
bash run_Snakebite-RNAseq.sh <module_or_rule> -np
```

The wrapper hardcodes `projectFolder`/`configFile` paths at the top and calls `snakemake` with
`--use-singularity`, a SLURM `--profile`, and `--restart-times 0 --retries 0` (retry/escalation
is instead handled per-rule, see "Resource escalation" below). This is a *project* wrapper -
copy it into each deployment and point it at that project's paths. It must be run from the
project folder (relative `config/` paths resolve against the current working directory, not
against `pipeline_folder`).

Modules (top-level rules in `rule all`): `reads`, `reference`, `preprocess`, `decontaminate`,
`align`, `assemble`, `quantify`, `report`. Each module can be run as a whole or as a single tool
within it using `<module>__<tool>__<action>` naming, e.g.:

```bash
bash run_Snakebite-RNAseq.sh align__star__map
bash run_Snakebite-RNAseq.sh assemble
```

Running a submodule automatically triggers its upstream dependencies.

## Architecture

### Config-file driven, not code-driven

Everything a user needs to change lives in `config/`, loaded once at the top of
`workflow/Snakefile`:

- `config/config.yaml` - paths to the other config files, `pipeline_folder`, temp-storage
  fallback chain, and `resource_sets` (named SLURM resource profiles: `small`, `medium`,
  `medium_nvme`, `large`, `large_nvme`, `highmem`, `longrun`, `highmem_longrun`, `full_node` -
  shared tier names with `Snakebite-Holoruminant-MetaG`/`Snakebite-miRNA`).
- `config/features.yaml` - the reference genome/annotation, and `decontamination:`, a
  name->path map of sequential decontamination references (empty/missing = decontamination
  skipped entirely).
- `config/params.yaml` - fastp/STAR/StringTie/featureCounts/RSEM parameters, including
  `align.star.mode` ("meta" for metatranscriptome studies - runs STAR without annotation, see
  below).
- `config/escalation.yaml` - per-rule resource tier, keyed by full rule name. Rules not listed
  default to `["small", "medium", "large"]` (`get_escalation_order` in
  `workflow/rules/helpers/__functions__.smk`) - in practice every rule is listed with a single
  best-fit tier (no rule currently uses real multi-tier retry escalation).
- `config/samples.tsv` - sample sheet: `sample_id`, `library_id` (one row per sequencing lane),
  `read1`, `read2`.
- `config/.docker.yml` - one container image per tool.
- `config/profiles/<cluster>/config.yaml` - SLURM executor profile; copy `Roihu/` for a new
  cluster. The default `resource_sets` are for CSC's Roihu, where `--gres=nvme` exists only on
  hugemem, so `nvme` stays empty and the local disk is `$TMPDIR`.

`workflow/Snakefile` parses `samples.tsv` into `SAMPLES` (unique `sample_id`s) and
`SAMPLE_LIBRARY` (`(sample_id, library_id)` pairs), and `features["decontamination"]` into
`CONTAM_NAMES` (order-preserving; empty list disables decontamination).

### Module layout under `workflow/rules/`

Each module is a directory with `__main__.smk` (includes the module's rule files and defines
the aggregate `rule <module>:`), `__functions__.smk` where the module has wildcard/lookup
helpers, and one `.smk` file per tool. `workflow/Snakefile` includes modules in pipeline order:
`helpers` -> `reads` -> `reference` -> `preprocess` -> `decontaminate` -> `align` -> `assemble`
-> `quantify` -> `report`.

`workflow/rules/folders.smk` centralizes every output-path constant (`TRIMMED`, `GENOME_STAR`,
`STRINGTIE`, etc., all `pathlib.Path` objects under `results/<module>/<tool>/...`) - check here
first when tracing where a rule's output goes. Every module gets its own `results/<module>/`
root, and a rule's `log:`/`benchmark:` files live *inside* that same tool folder (benchmarks
under a further `benchmark/` subfolder). The rendered report is the one exception, same as both
sibling pipelines: `Rreports/finalReport.html`, a separate, un-nested convention.

### The read-processing chain

Reads flow through the modules in a specific, mostly-linear chain:

```
raw lane (read1+read2) -> reads__fastqc__run / reads__multiqc__run (QC only, not in the chain)
raw lanes              -> preprocess__concatenate__run (lanes -> one FASTQ pair per sample)
                        -> preprocess__fastp__run (trim the concatenated pair)
                        -> decontaminate__star__map, once per reference in CONTAM_NAMES, chained
                           (each reference's unmapped-mate output feeds the next; skipped
                           entirely if CONTAM_NAMES is empty)
                        -> align__star__map (final genome+transcriptome alignment)
                        -> assemble__stringtie__run -> assemble__stringtie__merge
                        -> quantify__featurecounts__reference / __stringtie_merged
```

Note the trim order: concatenate *then* trim (the reverse of `Snakebite-miRNA`, which trims
each lane before concatenating) - this is this pipeline's own, deliberately-preserved order,
not a copy-paste of the sibling's.

`decontaminate/__functions__.smk`'s `get_final_read1`/`get_final_read2` resolve to either the
last reference's decontaminated output or straight to the trimmed (fastp) reads if
`CONTAM_NAMES` is empty - `align__star__map` calls these directly rather than depending on a
materialized "final" copy, since it's the only consumer.

### Metatranscriptome mode

`params["align"]["star"]["mode"] == "meta"` switches both `reference__star_index__genome` and
`align__star__map` to a Python-level `if`/`else` that defines a *different* rule body (same
pattern already used for optional decontamination) - the meta branch runs STAR without
`--sjdbGTFfile`/`--quantMode TranscriptomeSAM` and has no `transcriptome_bam` output, since STAR
never produces one in that mode. `quantify__rsem__calculate_expression` references that BAM
path directly (not via `rules.align__star__map.output.transcriptome_bam`) so that meta-mode
runs don't hit a parse-time `AttributeError` even when RSEM is never actually requested.

### Resource escalation (retry-with-bigger-resources)

Identical mechanism to the sibling pipelines: every rule uses `esc(key, rule_name)` /
`esc_val(key, rule_name)` (from `workflow/rules/helpers/__functions__.smk`) instead of
hardcoding `resources:`/`threads:`. `retries:` on each rule is
`len(get_escalation_order(rule_name))`, and a loop at the bottom of `workflow/Snakefile` forces
`r.retries = 0` for any rule whose escalation list has <=1 entry. Every STAR rule additionally
sets `--outTmpDir` to a rule/sample-local scratch path (cleared before each attempt) so retries
under SLURM's normal high-parallelism (`--jobs 300`) don't collide on or get blocked by a stale
temp dir.

### Temporary storage

No rule reads `tmp_storage`, `nvme_storage` or `shm_storage` at the moment, and `PROJECT_TMP` in
`folders.smk` isn't used either. STAR's `--outTmpDir` and featureCounts' temporary files are next
to their results, so no rule depends on the size of `$TMPDIR` (20 GiB in Roihu's `small`).

### RSEM stays opt-in

`reference__rsem_index__build` and `quantify__rsem__calculate_expression` are real, runnable
rules, callable directly (`bash run_Snakebite-RNAseq.sh reference__rsem_index__build`), but
excluded from the default `reference`/`quantify` module targets - same "available, not default"
precedent as `Snakebite-miRNA`'s bowtie-only genome-alignment branch.

### Report

`workflow/scripts/R/workflow-report.Rmd` is one large, pipeline-specific report (not split into
per-module Rmds like the sibling pipelines - it predates that convention and hasn't been broken
up). It expects `results/reads/multiqc_R1`, `results/preprocess/concatenate/multiqc_R1`,
`results/preprocess/fastp/multiqc_R1` (and their R2 counterparts), the featureCounts output
under `results/quantify/featurecounts/reference/`, and - if decontamination is active - a
`{sample_id}_R1.wcl` read-count file per reference under
`results/decontaminate/star/<reference>/` (written by `decontaminate__star__map`, consumed by
the report's per-reference decontamination-yield chart). It derives the configured
decontamination references itself by reading `config/features.yaml` as plain text and pulling
the keys under `decontamination:` (see the report rule's `features_file` param) - it does not
use a YAML parser, matching how the rest of the script already treats config as text rather
than adding a new R dependency.

### Version string

The pipeline version shown in logs is parsed out of the `CHANGELOG` file at Snakemake startup
(`get_version_from_changelog()` in `workflow/Snakefile`), matching a `Version X.Y.*:` header
followed by a leading patch number - bump `CHANGELOG`, not a separate version constant, when
releasing.
