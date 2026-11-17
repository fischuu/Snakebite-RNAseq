rule report__workflow__render:
    """Render the final HTML report (R Markdown)"""
    input:
        script=WORKFLOW_REPORT_R,
        mqc_raw=rules.reads__multiqc__run.output.r1,
        mqc_concatenated=rules.preprocess__multiqc__concatenated.output.r1,
        mqc_trimmed=rules.preprocess__multiqc__trimmed.output.r1,
        fc_reference=expand(str(QUANT_FEATURECOUNTS_REFERENCE / "{sample_id}_fc.txt"), sample_id=SAMPLES),
        fc_stringtie=expand(str(QUANT_FEATURECOUNTS_STRINGTIE / "{sample_id}_fc.txt"), sample_id=SAMPLES),
    output:
        PIPELINE_REPORT / "finalReport.html"
    log:
        PIPELINE_REPORT / "finalReport.log"
    benchmark:
        PIPELINE_REPORT / "benchmark" / "finalReport.tsv"
    params:
        project_folder=WD,
        pipeline_folder=config["pipeline_folder"],
        pipeline_config=config["params-file"],
        features_file=config["features-file"],
        samplesheet=config["sample-file"],
        sampleinfo=config["sampleinfo-file"],
        ref_annot_file=features["references"]["annotation"],
    threads: esc("cpus", "report__workflow__render")
    resources:
        runtime=esc("runtime", "report__workflow__render"),
        mem_mb=esc("mem_mb", "report__workflow__render"),
        cpus_per_task=esc("cpus", "report__workflow__render"),
        slurm_partition=esc("partition", "report__workflow__render"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'report__workflow__render')['nvme']}",
        attempt=get_attempt,
    container:
        docker["r_report"]
    shell:
        """
        exec > {log} 2>&1
        R -e "projFolder <- '{params.project_folder}'; \
              pipelineFolder <- '{params.pipeline_folder}'; \
              pipelineConfig <- '{params.pipeline_config}'; \
              features_file <- '{params.features_file}'; \
              refAnnot.file <- '{params.ref_annot_file}'; \
              samplesheet.file <- '{params.samplesheet}'; \
              sampleinfo.file <- '{params.sampleinfo}'; \
              snakemake <- TRUE; \
              rmarkdown::render('{input.script}', output_file='{output}')"
        """
