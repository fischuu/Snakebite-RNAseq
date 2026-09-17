include: "__functions__.smk"
include: "concatenate.smk"
include: "fastp.smk"
include: "fastqc.smk"

rule preprocess:
    """Concatenate lanes per sample, trim adapters, and QC the results"""
    input:
        expand(str(TRIMMED / "{sample_id}_R1.trimmed.fastq.gz"), sample_id=SAMPLES),
        expand(str(TRIMMED / "{sample_id}_R2.trimmed.fastq.gz"), sample_id=SAMPLES),
        rules.preprocess__multiqc__concatenated.output,
        rules.preprocess__multiqc__trimmed.output,
