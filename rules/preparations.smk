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