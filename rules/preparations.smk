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
    
rule cutadapt_trim_reads:
    """
    Trim adapters and low quality reads (CUTADAPT).
    """
    input:
        ["%s/%s/FASTQ/{samples}_1.fastq.gz" % (config["project-folder"], config["species"]),
         "%s/%s/FASTQ/{samples}_2.fastq.gz" % (config["project-folder"], config["species"])]
    output:
        first="%s/%s/cutadapt/{samples}_cutadapt_R1.fastq.gz" % (config["project-folder"], config["species"]),
        second="%s/%s/cutadapt/{samples}_cutadapt_R2.fastq.gz" % (config["project-folder"], config["species"])
    params:
        phread_score=config["params"]["cutadapt"]["phread_score"],
        adapter_file_R1=config["params"]["cutadapt"]["adapter_R1"],
        adapter_file_R2=config["params"]["cutadapt"]["adapter_R2"],
        pair_filter=config["params"]["cutadapt"]["pair_filter"],
        min_length=config["params"]["cutadapt"]["min_length"]
    log:
        "%s/%s/logs/cutadapt.{samples}.log" % (config["project-folder"], config["species"])
    benchmark:
        "%s/%s/benchmark/cutadapt.{samples}.benchmark.tsv" % (config["project-folder"], config["species"])
    threads: 1
    shell:"""
        output_dir=$(dirname {output.first})
        [ ! -d \"$output_dir\" ] && mkdir -p $output_dir

      	#******PARAMETERS*****
      	# -q : threshold used for quality trimming
      	# -a : path to file containing adapter sequence that might be ligated 3' end of the first read
      	# -A : path to file containing adapter sequence that might be ligated 3' end of the second read
      	# -o : path to output file for first read
      	# -p : path to output file for second read
      	# --minimum-length : reads shorter than this length are discarded
      	# --pair-filter : if "any" then the pair is discarded if one of the reads meet the filtering criterium (min length)

      	#******IN FILES*******
      	#gzipped or unzipped fastq forward and reverse reads

      	#*****OUT FILES*******

      	printf \"Removing low quality reads and Illumina adapter from the reads\\n\"

      	cutadapt -q {params.phread_score} -a file:{params.adapter_file_R1} -A file:{params.adapter_file_R2} -o {output.first} -p {output.second} --minimum-length={params.min_length} --pair-filter={params.pair_filter} {input} > {log} 2>&1    
	"""
