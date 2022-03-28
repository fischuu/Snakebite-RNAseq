# vim: set filetype=sh :
import pandas as pd
from snakemake.utils import validate, min_version

report: "report/workflow.rst"

##### RNASeq-snakemake pipeline #####
##### Daniel Fischer (daniel.fischer@luke.fi)
##### Natural Resources Institute Finland (Luke)
##### Version: 0.1.2

##### set minimum snakemake version #####
#min_version("6.0")

##### load config and sample sheets #####

rawsamples = pd.read_table(config["rawsamples"], header=None)[0].tolist()
samples = pd.read_table(config["samples"], header=None)[0].tolist()

workdir: config["project-folder"]

wildcard_constraints:
    rawsamples="|".join(rawsamples),
    samples="|".join(samples)

##### Functionm definitions #####

def make_fq_alignment_pairs(wildcards):
    sample_name = wildcards.samples
    if config["contamination-folder"] == "":
      r1= config["project-folder"]+"/FASTQ/TRIMMED/" + sample_name + "_R1.trimmed.fastq.gz"
      r2= config["project-folder"]+"/FASTQ/TRIMMED/" + sample_name + "_R2.trimmed.fastq.gz"
    else:
      r1= config["project-folder"]+"/FASTQ/DECONTAMINATED/" + sample_name + "_R1.decontaminated.fastq.gz"
      r2= config["project-folder"]+"/FASTQ/DECONTAMINATED/" + sample_name + "_R2.decontaminated.fastq.gz"
    result=[r1, r2]
    return result

##### run complete pipeline #####

rule all:
    input:
      # Step1-Preparations
        "%s/chrName.txt" % (config["star-index"]),
        "%s/Reference/RSEM.chrlist" % (config["project-folder"]),
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
#        expand(["%s/FASTQ/DECONTAMINATED/{samples}_R1.decontaminated.fastq.gz" % (config["project-folder"]),"%s/FASTQ/DECONTAMINATED/{samples}_R2.decontaminated.fastq.gz" % (config["project-folder"])], samples=samples),
      # Step4-Alignment
        expand("%s/BAM/{samples}.bam" % (config["project-folder"]), samples=samples),
      # Step5-TranscriptomeAssembly
        "%s/Stringtie/merged_STRG.gtf" % (config["project-folder"]),
      # Step6-Quantification
        expand("%s/RSEM/{samples}.genes.results" % (config["project-folder"]), samples=samples),
        expand("%s/FeatureCounts/Annot/{samples}_fc.txt" % (config["project-folder"]), samples=samples),
        expand("%s/FeatureCounts/Stringtie/{samples}_fc.txt" % (config["project-folder"]), samples=samples),
      # Step7-Reporting
        "%s/finalReport.html" % (config["project-folder"])
### setup report #####
report: "report/workflow.rst"

##### load rules #####
include: "rules/Step1-Preparations"
include: "rules/Step2-Preprocessing"
include: "rules/Step3-QC"
include: "rules/Step3b-Decontamination"
include: "rules/Step4-Alignment"
include: "rules/Step5-TranscriptomeAssembly"
include: "rules/Step6-Quantification"
include: "rules/Step7-Reporting"
