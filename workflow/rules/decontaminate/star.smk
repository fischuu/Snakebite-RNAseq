rule decontaminate__star_index__build:
    """Build the STAR index for one decontamination reference"""
    input:
        fasta=lambda wc: features["decontamination"][wc.contam],
    output:
        DECONTAM_INDEX / "{contam}" / "chrName.txt"
    log:
        DECONTAM_INDEX / "{contam}.log"
    benchmark:
        DECONTAM_INDEX / "benchmark" / "{contam}.tsv"
    params:
        index=lambda wc: str(DECONTAM_INDEX / wc.contam),
        tmpdir=lambda wc: str(DECONTAM_INDEX / wc.contam / "tmp"),
        limitGenomeGenerateRAM=params["align"]["star"]["limitGenomeGenerateRAM"],
    threads: esc("cpus", "decontaminate__star_index__build")
    resources:
        runtime=esc("runtime", "decontaminate__star_index__build"),
        mem_mb=esc("mem_mb", "decontaminate__star_index__build"),
        cpus_per_task=esc("cpus", "decontaminate__star_index__build"),
        slurm_partition=esc("partition", "decontaminate__star_index__build"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontaminate__star_index__build')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontaminate__star_index__build"))
    container:
        docker["star"]
    shell:
        """
        exec > {log} 2>&1
        rm -rf {params.tmpdir}
        STAR --runThreadN {threads} --outTmpDir {params.tmpdir} \
             --runMode genomeGenerate --genomeDir {params.index} \
             --genomeFastaFiles {input.fasta} \
             --limitGenomeGenerateRAM {params.limitGenomeGenerateRAM}
        """

rule decontaminate__star__map:
    """Map one sample against one decontamination reference and keep the unmapped pair
    (part of the sequential decontamination chain)"""
    input:
        index=DECONTAM_INDEX / "{contam}" / "chrName.txt",
        r1=get_decontam_chain_r1,
        r2=get_decontam_chain_r2,
    output:
        bam=DECONTAM_STAR / "{contam}" / "{sample_id}.bam",
        r1=DECONTAM_STAR / "{contam}" / "{sample_id}_R1.fastq.gz",
        r2=DECONTAM_STAR / "{contam}" / "{sample_id}_R2.fastq.gz",
        # Read count remaining after this reference, consumed by the report's
        # per-reference decontamination-yield chart (report__workflow__render).
        wcl=DECONTAM_STAR / "{contam}" / "{sample_id}_R1.wcl",
        logdir=directory(DECONTAM_STAR / "{contam}" / "{sample_id}"),
    log:
        DECONTAM_STAR / "{contam}" / "{sample_id}.log"
    benchmark:
        DECONTAM_STAR / "benchmark" / "{contam}" / "{sample_id}.tsv"
    params:
        index=lambda wc: str(DECONTAM_INDEX / wc.contam),
        tmpdir=lambda wc: str(DECONTAM_STAR / wc.contam / wc.sample_id / "tmp"),
        prefix=lambda wc: str(DECONTAM_STAR / wc.contam / wc.sample_id / f"{wc.sample_id}_"),
    threads: esc("cpus", "decontaminate__star__map")
    resources:
        runtime=esc("runtime", "decontaminate__star__map"),
        mem_mb=esc("mem_mb", "decontaminate__star__map"),
        cpus_per_task=esc("cpus", "decontaminate__star__map"),
        slurm_partition=esc("partition", "decontaminate__star__map"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontaminate__star__map')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontaminate__star__map"))
    container:
        docker["star"]
    shell:
        """
        exec > {log} 2>&1
        rm -rf {params.tmpdir}
        mkdir -p {output.logdir}

        STAR --genomeDir {params.index} --outTmpDir {params.tmpdir} \
             --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
             --runThreadN {threads} --outReadsUnmapped Fastx \
             --outFileNamePrefix {params.prefix}

        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
        wc -l < {params.prefix}Unmapped.out.mate1 > {output.wcl}
        gzip -c {params.prefix}Unmapped.out.mate1 > {output.r1}
        gzip -c {params.prefix}Unmapped.out.mate2 > {output.r2}
        """
