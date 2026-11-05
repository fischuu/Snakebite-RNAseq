# Snakebite-RNAseq

A Snakemake pipeline for paired-end bulk RNA-seq: per-lane FastQC, lane concatenation,
adapter/quality trimming (fastp), optional sequential decontamination against configurable
reference genomes (STAR), genome+transcriptome alignment (STAR, with a metatranscriptome mode),
transcript assembly and cross-sample merging (StringTie), gene/transcript quantification
(featureCounts; RSEM available on demand), and one consolidated HTML report. It runs on HPC
(SLURM) with every tool delivered as a Singularity/Apptainer container, and follows the same
module/escalation/docker architecture as
[Snakebite-Holoruminant-MetaG](https://github.com/fischuu/Snakebite-Holoruminant-MetaG) and
[Snakebite-miRNA](https://github.com/fischuu/Snakebite-miRNA).

# Installation

```
git clone git@github.com:fischuu/Snakebite-RNAseq.git
```

Or via HTTPS:

```
git clone https://github.com/fischuu/Snakebite-RNAseq.git
```

Test data is not part of this repository and needs to be downloaded separately.

# Running the pipeline

Copy `run_Snakebite-RNAseq.sh` into your project folder, point `projectFolder` at your project
and `pipeline_folder` (in `config/config.yaml`) at this checkout, then run it *from the project
folder*:

```bash
bash run_Snakebite-RNAseq.sh                  # full pipeline
bash run_Snakebite-RNAseq.sh align -np        # dry run of a single module
```

Modules (top-level rules in `rule all`): `reads`, `reference`, `preprocess`, `decontaminate`,
`align`, `assemble`, `quantify`, `report`. Each module can be run as a whole or as a single tool
within it using `<module>__<tool>__<action>` naming. Running a submodule automatically triggers
its upstream dependencies.

# Included tools

 * STAR
 * fastp
 * FastQC / MultiQC
 * StringTie
 * subread (featureCounts)
 * RSEM (available on demand, not part of the default `quantify` target)

# Configuration

 * `config/config.yaml` -- paths to the other config files and the SLURM resource tiers
   (`resource_sets`). The defaults are for CSC's Roihu. On another system, change them, and the
   limits in `tests/resources/check_resource_sets.py`, which checks them.
 * `config/features.yaml` -- the reference genome/annotation, and the (optional) sequential
   decontamination references.
 * `config/params.yaml` -- fastp/STAR/StringTie/featureCounts/RSEM parameters, including the
   metatranscriptome-mode toggle.
 * `config/escalation.yaml` -- per-rule resource tier.
 * `config/samples.tsv` -- sample sheet: `sample_id`, `library_id` (one row per sequencing
   lane), `read1`, `read2`.
 * `config/.docker.yml` -- one container image per tool.
 * `config/profiles/<cluster>/config.yaml` -- SLURM executor profile; copy `Roihu/` for a new
   cluster.

# Current status

Under development and possibly not yet stable for production.
