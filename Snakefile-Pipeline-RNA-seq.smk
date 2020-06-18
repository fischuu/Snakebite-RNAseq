# vim: set filetype=sh :
import pandas as pd
from snakemake.utils import validate, min_version

report: "report/workflow.rst"

##### set minimum snakemake version #####
#min_version("5.6")

##### load config and sample sheets #####

rawsamples = pd.read_table(config["rawsamples"], header=None)[0].tolist()
samples = pd.read_table(config["samples"], header=None)[0].tolist()

wildcard_constraints:
    rawsamples="|".join(rawsamples),
    samples="|".join(samples)

##### run complete pipeline #####

rule all:
    input:
      # Step1-Preparations
        "%s/chrName.txt" % (config["star-index"]),
        "%s/Reference/RSEM/" % (config["project-folder"]),
        expand("%s/FASTQ/CONCATENATED/{samples}_R1.fastq.gz" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/CONCATENATED/{samples}_R2.fastq.gz" % (config["project-folder"]), samples=samples),
      # Step2-Preprocessing
        "%s/chrName.txt" % (config["star-index"]),
        expand("%s/FASTQ/TRIMMED/{samples}_R1.trimmed.fastq.gz" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/{samples}_R2.trimmed.fastq.gz" % (config["project-folder"]), samples=samples),
      # Step3-QC
        "%s/QC/RAW/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/RAW/multiqc_R2/" % (config["project-folder"]),
        "%s/QC/CONCATENATED/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/CONCATENATED/multiqc_R2/" % (config["project-folder"]),
        "%s/QC/TRIMMED/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/TRIMMED/multiqc_R2/" % (config["project-folder"]),
      # Step4-Alignment
        expand("%s/BAM/{samples}.bam" % (config["project-folder"]), samples=samples),
      # Step5-Quantification
        expand("%s/RSEM/{samples}.genes.results" % (config["project-folder"]), samples=samples),
        expand("%s/FeatureCounts/{samples}_fc.txt" % (config["project-folder"]), samples=samples),

### setup report #####
report: "report/workflow.rst"

##### load rules #####
include: "rules/Step1-Preparations"
include: "rules/Step2-Preprocessing"
include: "rules/Step3-QC"
include: "rules/Step4-Alignment"
include: "rules/Step5-Quantification"
