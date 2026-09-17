rule reference__star_index__genome:
    """Build the STAR index for the reference genome"""
    input:
        fasta=features["references"]["genome"],
        annotation=features["references"]["annotation"],
    output:
        STAR_INDEX_GENOME / "chrName.txt"
    log:
        STAR_INDEX_GENOME / "star_index.log"
    benchmark:
        STAR_INDEX_GENOME / "benchmark.tsv"
    params:
        index=str(STAR_INDEX_GENOME),
        tmpdir=str(STAR_INDEX_GENOME / "tmp"),
        mode=params["align"]["star"]["mode"],
        overhang=params["align"]["star"]["overhang"],
        genomeSAindexNbases=params["align"]["star"]["genomeSAindexNbases"],
        genomeChrBinNbits=params["align"]["star"]["genomeChrBinNbits"],
        limitGenomeGenerateRAM=params["align"]["star"]["limitGenomeGenerateRAM"],
        sparseD=params["align"]["star"]["sparseD"],
    threads: esc("cpus", "reference__star_index__genome")
    resources:
        runtime=esc("runtime", "reference__star_index__genome"),
        mem_mb=esc("mem_mb", "reference__star_index__genome"),
        cpus_per_task=esc("cpus", "reference__star_index__genome"),
        slurm_partition=esc("partition", "reference__star_index__genome"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__star_index__genome')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__star_index__genome"))
    container:
        docker["star"]
    shell:
        """
        exec > {log} 2>&1
        rm -rf {params.tmpdir}

        if [ "{params.mode}" = "meta" ]; then
            STAR --runThreadN {threads} --outTmpDir {params.tmpdir} \
                 --runMode genomeGenerate --genomeDir {params.index} \
                 --genomeFastaFiles {input.fasta} \
                 --limitGenomeGenerateRAM {params.limitGenomeGenerateRAM} \
                 --genomeSAindexNbases {params.genomeSAindexNbases} \
                 --genomeChrBinNbits {params.genomeChrBinNbits} \
                 --genomeSAsparseD {params.sparseD}
        else
            STAR --runThreadN {threads} --outTmpDir {params.tmpdir} \
                 --runMode genomeGenerate --genomeDir {params.index} \
                 --genomeFastaFiles {input.fasta} \
                 --sjdbGTFfile {input.annotation} \
                 --limitGenomeGenerateRAM {params.limitGenomeGenerateRAM} \
                 --genomeSAindexNbases {params.genomeSAindexNbases} \
                 --genomeChrBinNbits {params.genomeChrBinNbits} \
                 --genomeSAsparseD {params.sparseD} \
                 --sjdbOverhang {params.overhang}
        fi
        """
