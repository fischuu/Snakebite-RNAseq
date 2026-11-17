# RSEM is kept available for on-demand use but is not part of the default `quantify` target
# (its featureCounts-based siblings are); see reference/rsem_index.smk.
#
# Only meaningful in non-meta mode: align__star__map's meta-mode branch runs STAR without
# --quantMode TranscriptomeSAM (no annotation, no transcriptome), so this input path is
# referenced directly (not via rules.align__star__map.output) to avoid a parse-time
# AttributeError against that branch's output.

rule quantify__rsem__calculate_expression:
    """Calculate one sample's transcript/gene expression from its transcriptome alignment (RSEM)"""
    input:
        bam=GENOME_STAR / "{sample_id}" / "{sample_id}_Aligned.toTranscriptome.out.bam",
        index=rules.reference__rsem_index__build.output,
    output:
        QUANT_RSEM / "{sample_id}.genes.results"
    log:
        QUANT_RSEM / "{sample_id}.log"
    benchmark:
        QUANT_RSEM / "benchmark" / "{sample_id}.tsv"
    params:
        rsem_index=str(RSEM_INDEX / "RSEM"),
        out_prefix=lambda wc: str(QUANT_RSEM / wc.sample_id),
        paired=params["quantify"]["rsem"]["paired"],
    threads: esc("cpus", "quantify__rsem__calculate_expression")
    resources:
        runtime=esc("runtime", "quantify__rsem__calculate_expression"),
        mem_mb=esc("mem_mb", "quantify__rsem__calculate_expression"),
        cpus_per_task=esc("cpus", "quantify__rsem__calculate_expression"),
        slurm_partition=esc("partition", "quantify__rsem__calculate_expression"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'quantify__rsem__calculate_expression')['nvme']}",
        attempt=get_attempt,
    container:
        docker["rsem"]
    shell:
        """
        exec > {log} 2>&1
        rsem-calculate-expression --alignments {input.bam} -p {threads} {params.paired} \
                                   --fragment-length-mean 75 \
                                   {params.rsem_index} {params.out_prefix}
        """
