# Tests

Checks to run when you change the pipeline. They need no data, databases or real tools.

| Test | What it proves | How long |
|------|----------------|----------|
| `escalation/run_test.sh` | A failed job escalates to a bigger SLURM job, with one attempt per resource set and no retry inside a job. Needs Docker | About 5 minutes |
| `resources/check_resource_sets.py` | The resource sets of `config/config.yaml` fit the partitions of CSC's Roihu, and every set in `config/escalation.yaml` exists. Needs only Python with PyYAML | A second |

The escalation test uses the one-node SLURM cluster in `minicluster/`, the same as in Snakebite-MetaFungi and in the book
of Snakebite-Holoruminant-MetaG. `IMAGE=<image> bash tests/escalation/run_test.sh` uses an image you built already, for
example with another version of Snakemake.
