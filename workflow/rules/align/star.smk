if params["align"]["star"]["mode"] == "meta":

    rule align__star__map:
        """Map one sample's fully preprocessed reads to the reference genome (STAR, metatranscriptome mode: no annotation)"""
        input:
            index=rules.reference__star_index__genome.output,
            r1=get_final_read1,
            r2=get_final_read2,
        output:
            bam=GENOME_STAR / "{sample_id}.bam",
            logdir=directory(GENOME_STAR / "{sample_id}"),
        log:
            GENOME_STAR / "{sample_id}.log"
        benchmark:
            GENOME_STAR / "benchmark" / "{sample_id}.tsv"
        params:
            index=str(STAR_INDEX_GENOME),
            tmpdir=lambda wc: str(GENOME_STAR / wc.sample_id / "tmp"),
            prefix=lambda wc: str(GENOME_STAR / wc.sample_id / f"{wc.sample_id}_"),
        threads: esc("cpus", "align__star__map")
        resources:
            runtime=esc("runtime", "align__star__map"),
            mem_mb=esc("mem_mb", "align__star__map"),
            cpus_per_task=esc("cpus", "align__star__map"),
            slurm_partition=esc("partition", "align__star__map"),
            gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__star__map')['nvme']}",
            attempt=get_attempt,
        retries: len(get_escalation_order("align__star__map"))
        container:
            docker["star"]
        shell:
            """
            exec > {log} 2>&1
            rm -rf {params.tmpdir}
            mkdir -p {output.logdir}

            STAR --genomeDir {params.index} --outTmpDir {params.tmpdir} \
                 --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat \
                 --outSAMunmapped Within --outSAMtype BAM Unsorted SortedByCoordinate \
                 --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical \
                 --runThreadN {threads} --alignSoftClipAtReferenceEnds No \
                 --outFileNamePrefix {params.prefix}

            mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
            """

else:

    rule align__star__map:
        """Map one sample's fully preprocessed reads to the reference genome and transcriptome (STAR)"""
        input:
            index=rules.reference__star_index__genome.output,
            annotation=features["references"]["annotation"],
            r1=get_final_read1,
            r2=get_final_read2,
        output:
            bam=GENOME_STAR / "{sample_id}.bam",
            transcriptome_bam=GENOME_STAR / "{sample_id}" / "{sample_id}_Aligned.toTranscriptome.out.bam",
            logdir=directory(GENOME_STAR / "{sample_id}"),
        log:
            GENOME_STAR / "{sample_id}.log"
        benchmark:
            GENOME_STAR / "benchmark" / "{sample_id}.tsv"
        params:
            index=str(STAR_INDEX_GENOME),
            tmpdir=lambda wc: str(GENOME_STAR / wc.sample_id / "tmp"),
            prefix=lambda wc: str(GENOME_STAR / wc.sample_id / f"{wc.sample_id}_"),
        threads: esc("cpus", "align__star__map")
        resources:
            runtime=esc("runtime", "align__star__map"),
            mem_mb=esc("mem_mb", "align__star__map"),
            cpus_per_task=esc("cpus", "align__star__map"),
            slurm_partition=esc("partition", "align__star__map"),
            gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__star__map')['nvme']}",
            attempt=get_attempt,
        retries: len(get_escalation_order("align__star__map"))
        container:
            docker["star"]
        shell:
            """
            exec > {log} 2>&1
            rm -rf {params.tmpdir}
            mkdir -p {output.logdir}

            STAR --genomeDir {params.index} --outTmpDir {params.tmpdir} \
                 --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat \
                 --outFilterType BySJout --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
                 --outSAMstrandField intronMotif --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical \
                 --runThreadN {threads} --quantMode TranscriptomeSAM --outWigType bedGraph --outWigStrand Stranded \
                 --alignSoftClipAtReferenceEnds No --outFileNamePrefix {params.prefix}

            mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
            """
