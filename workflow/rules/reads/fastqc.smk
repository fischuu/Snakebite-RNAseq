rule reads__fastqc__run:
    """Quality control of one raw lane's paired FASTQ (FastQC)"""
    input:
        r1=get_raw_read1,
        r2=get_raw_read2,
    output:
        r1=READS / "{sample_id}.{library_id}_R1_fastqc.zip",
        r1_html=READS / "{sample_id}.{library_id}_R1_fastqc.html",
        r2=READS / "{sample_id}.{library_id}_R2_fastqc.zip",
        r2_html=READS / "{sample_id}.{library_id}_R2_fastqc.html",
    log:
        READS / "{sample_id}.{library_id}_fastqc.log"
    benchmark:
        READS / "benchmark/{sample_id}.{library_id}.tsv"
    params:
        outfolder=str(READS),
    threads: esc("cpus", "reads__fastqc__run")
    resources:
        runtime=esc("runtime", "reads__fastqc__run"),
        mem_mb=esc("mem_mb", "reads__fastqc__run"),
        cpus_per_task=esc("cpus", "reads__fastqc__run"),
        slurm_partition=esc("partition", "reads__fastqc__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reads__fastqc__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reads__fastqc__run"))
    container:
        docker["fastqc"]
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outfolder}

        fastqc -t {threads} -o {params.outfolder} --extract {input.r1}
        base=$(basename {input.r1})
        base=${{base%.fastq.gz}}
        base=${{base%.fastq}}
        mv {params.outfolder}/${{base}}_fastqc.zip {output.r1}
        mv {params.outfolder}/${{base}}_fastqc.html {output.r1_html}

        fastqc -t {threads} -o {params.outfolder} --extract {input.r2}
        base=$(basename {input.r2})
        base=${{base%.fastq.gz}}
        base=${{base%.fastq}}
        mv {params.outfolder}/${{base}}_fastqc.zip {output.r2}
        mv {params.outfolder}/${{base}}_fastqc.html {output.r2_html}
        """

rule reads__multiqc__run:
    """Aggregate FastQC reports for every raw lane, separately for R1 and R2 (MultiQC)"""
    input:
        r1=[READS / f"{sample_id}.{library_id}_R1_fastqc.zip" for sample_id, library_id in SAMPLE_LIBRARY],
        r2=[READS / f"{sample_id}.{library_id}_R2_fastqc.zip" for sample_id, library_id in SAMPLE_LIBRARY],
    output:
        r1=directory(READS / "multiqc_R1"),
        r2=directory(READS / "multiqc_R2"),
    log:
        r1=READS / "multiqc_R1.log",
        r2=READS / "multiqc_R2.log",
    benchmark:
        READS / "benchmark/multiqc.tsv"
    threads: esc("cpus", "reads__multiqc__run")
    resources:
        runtime=esc("runtime", "reads__multiqc__run"),
        mem_mb=esc("mem_mb", "reads__multiqc__run"),
        cpus_per_task=esc("cpus", "reads__multiqc__run"),
        slurm_partition=esc("partition", "reads__multiqc__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reads__multiqc__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reads__multiqc__run"))
    container:
        docker["multiqc"]
    shell:
        """
        multiqc -f --interactive -o {output.r1} {input.r1} &> {log.r1}
        multiqc -f --interactive -o {output.r2} {input.r2} &> {log.r2}
        """
