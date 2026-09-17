rule quantify__featurecounts__reference:
    """Quantify one sample's aligned reads against the reference annotation (featureCounts)"""
    input:
        bam=rules.align__star__map.output.bam,
        gtf=features["references"]["annotation"],
    output:
        QUANT_FEATURECOUNTS_REFERENCE / "{sample_id}_fc.txt"
    log:
        QUANT_FEATURECOUNTS_REFERENCE / "{sample_id}.log"
    benchmark:
        QUANT_FEATURECOUNTS_REFERENCE / "benchmark" / "{sample_id}.tsv"
    params:
        geneid=params["quantify"]["featurecounts"]["attributetype"],
        type=params["quantify"]["featurecounts"]["featuretype"],
        strand=params["quantify"]["featurecounts"]["strandtype"],
    threads: esc("cpus", "quantify__featurecounts__reference")
    resources:
        runtime=esc("runtime", "quantify__featurecounts__reference"),
        mem_mb=esc("mem_mb", "quantify__featurecounts__reference"),
        cpus_per_task=esc("cpus", "quantify__featurecounts__reference"),
        slurm_partition=esc("partition", "quantify__featurecounts__reference"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'quantify__featurecounts__reference')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("quantify__featurecounts__reference"))
    container:
        docker["subread"]
    shell:
        """
        exec > {log} 2>&1
        featureCounts -p -T {threads} \
                      -a {input.gtf} -o {output} \
                      -g {params.geneid} -t {params.type} -s {params.strand} \
                      {input.bam}
        """

rule quantify__featurecounts__stringtie_merged:
    """Quantify one sample's aligned reads against the merged StringTie transcript assembly (featureCounts)"""
    input:
        bam=rules.align__star__map.output.bam,
        gtf=rules.assemble__stringtie__merge.output,
    output:
        QUANT_FEATURECOUNTS_STRINGTIE / "{sample_id}_fc.txt"
    log:
        QUANT_FEATURECOUNTS_STRINGTIE / "{sample_id}.log"
    benchmark:
        QUANT_FEATURECOUNTS_STRINGTIE / "benchmark" / "{sample_id}.tsv"
    params:
        # Always transcript_id/transcript for the merged-assembly count, regardless of the
        # reference-annotation attribute/feature type above - StringTie's merged GTF only
        # carries transcript-level records.
        strand=params["quantify"]["featurecounts"]["strandtype"],
    threads: esc("cpus", "quantify__featurecounts__stringtie_merged")
    resources:
        runtime=esc("runtime", "quantify__featurecounts__stringtie_merged"),
        mem_mb=esc("mem_mb", "quantify__featurecounts__stringtie_merged"),
        cpus_per_task=esc("cpus", "quantify__featurecounts__stringtie_merged"),
        slurm_partition=esc("partition", "quantify__featurecounts__stringtie_merged"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'quantify__featurecounts__stringtie_merged')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("quantify__featurecounts__stringtie_merged"))
    container:
        docker["subread"]
    shell:
        """
        exec > {log} 2>&1
        featureCounts -p -T {threads} \
                      -a {input.gtf} -o {output} \
                      -g transcript_id -t transcript -s {params.strand} \
                      {input.bam}
        """
