include: "star.smk"

rule align:
    """Map every sample's fully preprocessed reads to the reference genome (STAR)"""
    input:
        expand(str(GENOME_STAR / "{sample_id}.bam"), sample_id=SAMPLES),
