# vim: set filetype=sh :
import pandas as pd
from snakemake.utils import validate, min_version
from multiprocessing import cpu_count

report: "report/workflow.rst"

##### set minimum snakemake version #####
#min_version("5.1.2")

##### load config and sample sheets #####

samples = pd.read_table(config["samples"], header=None)[0].tolist()

##### run complete pipeline #####

rule all:
    input:
      # Prepare the genome index for STAR
        "%s/chrName.txt" % (config["index"]),
      # Index the genome with samtools
        "%s.fai" % (config["genome"]),
      # Quality check raw data:
        "%s/QC/raw/multiqc/" % (config["project-folder"]),
      # Trim low quality reads
        expand("%s/FASTQ/cutadapt/{samples}_cutadapt.fastq.gz" % (config["project-folder"]), samples=samples),
      # Quality check for trimmed data
        "%s/QC/cutadapt/multiqc/" % (config["project-folder"]),
      # Map the reads
        expand("%s/BAM/{samples}.bam" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
      # Variant calling
        "%s/VCF/samtools_mpileup.bcf" % (config["project-folder"]),
        "%s/VCF/samtools_var.bcf" % (config["project-folder"]),
        "%s/VCF/samtools_variants.vcf" % (config["project-folder"]),
      # Quantify the reads
        expand("%s/GTF/FeatureCounts/{samples}_fc.txt" % (config["project-folder"]), samples=samples),
      # Create RSEM index
        "%s/Reference/RSEM/" % (config["project-folder"]),
      # Quantify with RSEM
        expand("%s/RSEM/{samples}.genes.results" % (config["project-folder"]), samples=samples)

### setup report #####
report: "report/workflow.rst"

##### load rules #####
include: "rules/preparations.smk"
#include: "rules/qc.smk"
#include: "rules/mapping.smk"
#inlcude: "rules/quantification.smk"


#include: "rules/star_create_index.smk"
include: "rules/fastqc_quality_control.smk"
include: "rules/multiqc_quality_control.smk"
include: "rules/cutadapt_trim_reads.smk"
include: "rules/fastqc_quality_control_cutadapt.smk"
include: "rules/multiqc_quality_control_cutadapt.smk"
include: "rules/star_map_reads.smk"
include: "rules/featureCounts_quantify.smk"
include: "rules/rsem_create_index.smk"
include: "rules/rsem_calculate_expressions.smk"
include: "rules/samtools_bam_tosortedBam.smk"
include: "rules/samtools_indexFasta.smk"
include: "rules/samtools_mpileup.smk"
include: "rules/samtools_callvariants.smk"
include: "rules/samtools_filtervariants.smk"