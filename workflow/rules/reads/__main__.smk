include: "__functions__.smk"
include: "fastqc.smk"

rule reads:
    """Run FastQC/MultiQC on every raw lane"""
    input:
        rules.reads__multiqc__run.output,
