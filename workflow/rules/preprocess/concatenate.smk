rule preprocess__concatenate__run:
    """Concatenate every lane belonging to one sample into a single paired FASTQ"""
    input:
        r1=get_concatenation_inputs_read1,
        r2=get_concatenation_inputs_read2,
    output:
        r1=temp(CONCATENATED / "{sample_id}_R1.fastq.gz"),
        r2=temp(CONCATENATED / "{sample_id}_R2.fastq.gz"),
    log:
        CONCATENATED / "{sample_id}.log"
    benchmark:
        CONCATENATED / "benchmark/{sample_id}.tsv"
    threads: esc("cpus", "preprocess__concatenate__run")
    resources:
        runtime=esc("runtime", "preprocess__concatenate__run"),
        mem_mb=esc("mem_mb", "preprocess__concatenate__run"),
        cpus_per_task=esc("cpus", "preprocess__concatenate__run"),
        slurm_partition=esc("partition", "preprocess__concatenate__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preprocess__concatenate__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__concatenate__run"))
    shell:
        """
        exec > {log} 2>&1
        echo "Concatenating read1: {input.r1}"
        cat {input.r1} > {output.r1}
        echo "Concatenating read2: {input.r2}"
        cat {input.r2} > {output.r2}
        """
