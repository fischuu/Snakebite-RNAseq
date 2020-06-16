# vim: set filetype=sh :

rule star_create_index:
    """
    Create genome index for reference genome (STAR).
    """
    input:
        index="%s" % (config["star-index"]),
        fasta="%s" % (config["genome"]),
        annot="%s" % (config["annot"])
    output:
        "%s/chrName.txt" % (config["star-index"])
    log:
        "%s/logs/STAR/createIndex.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/STAR/createIndex.benchmark.tsv" % (config["project-folder"])
    conda: "../envs/star.yaml"
    threads: lambda cores: cpu_count()
    shell:"""
        echo "Number of threads used:" {threads}
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {input.index} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.annot} --sjdbOverhang 49 2> {log}
    """
    
rule Step0_Concatenate_lanes:
    """
    Concatenate the demultiplexed fastq files (BASH).
    """
    input:
      R1=expand("%s/FASTQ/RAW/{rawsamples}_R1_001.fastq.gz" % (config["project-folder"]), rawsamples=rawsamples),
      R2=expand("%s/FASTQ/RAW/{rawsamples}_R2_001.fastq.gz" % (config["project-folder"]), rawsamples=rawsamples)
    output:
      R1="%s/FASTQ/CONCATENATED/{samples}_R1.merged.fastq.gz" % (config["project-folder"]),
      R2="%s/FASTQ/CONCATENATED/{samples}_R2.merged.fastq.gz" % (config["project-folder"])
    log:
        "%s/logs/Concatenate/catFastq_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Concatenate/{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 1
    params:
       infolder="%s/FASTQ/RAW" % (config["project-folder"]),
       outfolder="%s/FASTQ/CONCATENATED" % (config["project-folder"])
    shell:"""
        echo "Number of threads used:" {threads}
        mkdir -p {params.outfolder}
        cat {params.infolder}/{wildcards.samples}*_R1_001.fastq.gz > {output.R1} 2> {log}
        cat {params.infolder}/{wildcards.samples}*_R2_001.fastq.gz > {output.R2} 2> {log}
  	"""    
    
