rule preprocess__fastqc__concatenated:
    """Quality control of one sample's concatenated FASTQ (FastQC)"""
    input:
        r1=rules.preprocess__concatenate__run.output.r1,
        r2=rules.preprocess__concatenate__run.output.r2,
    output:
        r1=CONCATENATED / "{sample_id}_R1_fastqc.zip",
        r2=CONCATENATED / "{sample_id}_R2_fastqc.zip",
    log:
        CONCATENATED / "{sample_id}_fastqc.log"
    benchmark:
        CONCATENATED / "benchmark/fastqc.{sample_id}.tsv"
    params:
        outfolder=str(CONCATENATED),
    threads: esc("cpus", "preprocess__fastqc__concatenated")
    resources:
        runtime=esc("runtime", "preprocess__fastqc__concatenated"),
        mem_mb=esc("mem_mb", "preprocess__fastqc__concatenated"),
        cpus_per_task=esc("cpus", "preprocess__fastqc__concatenated"),
        slurm_partition=esc("partition", "preprocess__fastqc__concatenated"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preprocess__fastqc__concatenated')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__fastqc__concatenated"))
    container:
        docker["fastqc"]
    shell:
        """
        exec > {log} 2>&1
        fastqc -t {threads} -o {params.outfolder} --extract {input.r1}
        fastqc -t {threads} -o {params.outfolder} --extract {input.r2}
        """

rule preprocess__multiqc__concatenated:
    """Aggregate FastQC reports for every sample's concatenated FASTQ, separately for R1 and R2 (MultiQC)"""
    input:
        r1=expand(str(CONCATENATED / "{sample_id}_R1_fastqc.zip"), sample_id=SAMPLES),
        r2=expand(str(CONCATENATED / "{sample_id}_R2_fastqc.zip"), sample_id=SAMPLES),
    output:
        r1=directory(CONCATENATED / "multiqc_R1"),
        r2=directory(CONCATENATED / "multiqc_R2"),
    log:
        r1=CONCATENATED / "multiqc_R1.log",
        r2=CONCATENATED / "multiqc_R2.log",
    benchmark:
        CONCATENATED / "benchmark/multiqc.tsv"
    threads: esc("cpus", "preprocess__multiqc__concatenated")
    resources:
        runtime=esc("runtime", "preprocess__multiqc__concatenated"),
        mem_mb=esc("mem_mb", "preprocess__multiqc__concatenated"),
        cpus_per_task=esc("cpus", "preprocess__multiqc__concatenated"),
        slurm_partition=esc("partition", "preprocess__multiqc__concatenated"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preprocess__multiqc__concatenated')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__multiqc__concatenated"))
    container:
        docker["multiqc"]
    shell:
        """
        multiqc -f --interactive -o {output.r1} {input.r1} &> {log.r1}
        multiqc -f --interactive -o {output.r2} {input.r2} &> {log.r2}
        """

rule preprocess__fastqc__trimmed:
    """Quality control of one sample's trimmed FASTQ (FastQC)"""
    input:
        r1=rules.preprocess__fastp__run.output.r1,
        r2=rules.preprocess__fastp__run.output.r2,
    output:
        r1=TRIMMED / "{sample_id}_R1.trimmed_fastqc.zip",
        r2=TRIMMED / "{sample_id}_R2.trimmed_fastqc.zip",
    log:
        TRIMMED / "{sample_id}_fastqc.log"
    benchmark:
        TRIMMED / "benchmark/fastqc.{sample_id}.tsv"
    params:
        outfolder=str(TRIMMED),
    threads: esc("cpus", "preprocess__fastqc__trimmed")
    resources:
        runtime=esc("runtime", "preprocess__fastqc__trimmed"),
        mem_mb=esc("mem_mb", "preprocess__fastqc__trimmed"),
        cpus_per_task=esc("cpus", "preprocess__fastqc__trimmed"),
        slurm_partition=esc("partition", "preprocess__fastqc__trimmed"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preprocess__fastqc__trimmed')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__fastqc__trimmed"))
    container:
        docker["fastqc"]
    shell:
        """
        exec > {log} 2>&1
        fastqc -t {threads} -o {params.outfolder} --extract {input.r1}
        fastqc -t {threads} -o {params.outfolder} --extract {input.r2}
        """

rule preprocess__multiqc__trimmed:
    """Aggregate FastQC reports for every sample's trimmed FASTQ, separately for R1 and R2 (MultiQC)"""
    input:
        r1=expand(str(TRIMMED / "{sample_id}_R1.trimmed_fastqc.zip"), sample_id=SAMPLES),
        r2=expand(str(TRIMMED / "{sample_id}_R2.trimmed_fastqc.zip"), sample_id=SAMPLES),
    output:
        r1=directory(TRIMMED / "multiqc_R1"),
        r2=directory(TRIMMED / "multiqc_R2"),
    log:
        r1=TRIMMED / "multiqc_R1.log",
        r2=TRIMMED / "multiqc_R2.log",
    benchmark:
        TRIMMED / "benchmark/multiqc.tsv"
    threads: esc("cpus", "preprocess__multiqc__trimmed")
    resources:
        runtime=esc("runtime", "preprocess__multiqc__trimmed"),
        mem_mb=esc("mem_mb", "preprocess__multiqc__trimmed"),
        cpus_per_task=esc("cpus", "preprocess__multiqc__trimmed"),
        slurm_partition=esc("partition", "preprocess__multiqc__trimmed"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preprocess__multiqc__trimmed')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__multiqc__trimmed"))
    container:
        docker["multiqc"]
    shell:
        """
        multiqc -f --interactive -o {output.r1} {input.r1} &> {log.r1}
        multiqc -f --interactive -o {output.r2} {input.r2} &> {log.r2}
        """
