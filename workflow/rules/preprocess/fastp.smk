rule preprocess__fastp__run:
    """Trim adapters and low-quality bases from one sample's concatenated FASTQ (fastp)"""
    input:
        r1=rules.preprocess__concatenate__run.output.r1,
        r2=rules.preprocess__concatenate__run.output.r2,
    output:
        r1=TRIMMED / "{sample_id}_R1.trimmed.fastq.gz",
        r2=TRIMMED / "{sample_id}_R2.trimmed.fastq.gz",
    log:
        TRIMMED / "{sample_id}_fastp.log"
    benchmark:
        TRIMMED / "benchmark/{sample_id}.tsv"
    params:
        trim_front1=params["preprocess"]["fastp"]["trim_front1"],
        trim_tail1=params["preprocess"]["fastp"]["trim_tail1"],
        trim_front2=params["preprocess"]["fastp"]["trim_front2"],
        trim_tail2=params["preprocess"]["fastp"]["trim_tail2"],
    threads: esc("cpus", "preprocess__fastp__run")
    resources:
        runtime=esc("runtime", "preprocess__fastp__run"),
        mem_mb=esc("mem_mb", "preprocess__fastp__run"),
        cpus_per_task=esc("cpus", "preprocess__fastp__run"),
        slurm_partition=esc("partition", "preprocess__fastp__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preprocess__fastp__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__fastp__run"))
    container:
        docker["fastp"]
    shell:
        """
        exec > {log} 2>&1
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
              --trim_front1 {params.trim_front1} --trim_tail1 {params.trim_tail1} \
              --trim_front2 {params.trim_front2} --trim_tail2 {params.trim_tail2} \
              --thread {threads}
        """
