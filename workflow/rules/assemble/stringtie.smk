rule assemble__stringtie__run:
    """Assemble one sample's transcripts from its genome alignment (StringTie)"""
    input:
        bam=rules.align__star__map.output.bam,
        annotation=features["references"]["annotation"],
    output:
        STRINGTIE / "{sample_id}.stringtie.gtf"
    log:
        STRINGTIE / "{sample_id}.log"
    benchmark:
        STRINGTIE / "benchmark" / "{sample_id}.tsv"
    threads: esc("cpus", "assemble__stringtie__run")
    resources:
        runtime=esc("runtime", "assemble__stringtie__run"),
        mem_mb=esc("mem_mb", "assemble__stringtie__run"),
        cpus_per_task=esc("cpus", "assemble__stringtie__run"),
        slurm_partition=esc("partition", "assemble__stringtie__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'assemble__stringtie__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__stringtie__run"))
    container:
        docker["stringtie"]
    shell:
        """
        exec > {log} 2>&1
        stringtie {input.bam} -p {threads} -G {input.annotation} -v -o {output}
        """

rule assemble__stringtie__compose_merge:
    """List every sample's assembled GTF for stringtie --merge"""
    input:
        files=expand(str(STRINGTIE / "{sample_id}.stringtie.gtf"), sample_id=SAMPLES)
    output:
        txt=STRINGTIE / "stringtie_gtfs.txt"
    log:
        STRINGTIE / "compose_merge.log"
    benchmark:
        STRINGTIE / "benchmark" / "compose_merge.tsv"
    params:
        folder=str(STRINGTIE),
        script=os.path.join(SCRIPT_FOLDER, "composeMerge.sh"),
    threads: esc("cpus", "assemble__stringtie__compose_merge")
    resources:
        runtime=esc("runtime", "assemble__stringtie__compose_merge"),
        mem_mb=esc("mem_mb", "assemble__stringtie__compose_merge"),
        cpus_per_task=esc("cpus", "assemble__stringtie__compose_merge"),
        slurm_partition=esc("partition", "assemble__stringtie__compose_merge"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'assemble__stringtie__compose_merge')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__stringtie__compose_merge"))
    shell:
        """
        exec > {log} 2>&1
        {params.script} {params.folder} {output.txt}
        """

rule assemble__stringtie__merge:
    """Merge every sample's assembled transcripts into one consensus GTF (StringTie --merge)"""
    input:
        gtfs=rules.assemble__stringtie__compose_merge.output.txt,
        annotation=features["references"]["annotation"],
    output:
        STRINGTIE / "merged_STRG.gtf"
    log:
        STRINGTIE / "merge.log"
    benchmark:
        STRINGTIE / "benchmark" / "merge.tsv"
    params:
        tpm=params["assemble"]["stringtie"]["tpm"],
    threads: esc("cpus", "assemble__stringtie__merge")
    resources:
        runtime=esc("runtime", "assemble__stringtie__merge"),
        mem_mb=esc("mem_mb", "assemble__stringtie__merge"),
        cpus_per_task=esc("cpus", "assemble__stringtie__merge"),
        slurm_partition=esc("partition", "assemble__stringtie__merge"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'assemble__stringtie__merge')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("assemble__stringtie__merge"))
    container:
        docker["stringtie"]
    shell:
        """
        exec > {log} 2>&1
        stringtie --merge -G {input.annotation} -F 0 -T {params.tpm} -o {output} {input.gtfs}
        """
