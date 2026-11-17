# RSEM is kept available for on-demand use (see quantify/rsem.smk) but is not part of the
# default `reference`/`quantify` targets.

rule reference__rsem_index__build:
    """Build the RSEM reference index"""
    input:
        fasta=features["references"]["genome"],
        annotation=features["references"]["annotation"],
    output:
        RSEM_INDEX / "RSEM.chrlist"
    log:
        RSEM_INDEX / "rsem_index.log"
    benchmark:
        RSEM_INDEX / "benchmark.tsv"
    params:
        out=str(RSEM_INDEX / "RSEM"),
    threads: esc("cpus", "reference__rsem_index__build")
    resources:
        runtime=esc("runtime", "reference__rsem_index__build"),
        mem_mb=esc("mem_mb", "reference__rsem_index__build"),
        cpus_per_task=esc("cpus", "reference__rsem_index__build"),
        slurm_partition=esc("partition", "reference__rsem_index__build"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__rsem_index__build')['nvme']}",
        attempt=get_attempt,
    container:
        docker["rsem"]
    shell:
        """
        exec > {log} 2>&1
        rsem-prepare-reference -p {threads} --gtf {input.annotation} --star {input.fasta} {params.out}
        """
