import pandas as pd
from snakemake.utils import validate, min_version
from multiprocessing import cpu_count
import glob
import re
import os
import sys
import yaml

report: "report/workflow.rst"

##### RNASeq-snakemake pipeline #####
##### Daniel Fischer (daniel.fischer@luke.fi)
##### Natural Resources Institute Finland (Luke)
##### Version: 0.2
version = "0.2"

##### set minimum snakemake version #####
#min_version("6.0")

##### load config and sample sheets #####

workdir: config["project-folder"]

##### Fill the configuration lines for relative paths
if config["project-folder"][-1] == '/':
   config["project-folder"]=config["project-folder"][:-1]
   
if(config["pipeline-config"][0]!='/'):
    config["pipeline-config"] = config["project-folder"] + '/' + config["pipeline-config"]

if(config["server-config"][0]!='/'):
    config["server-config"] = config["project-folder"] + '/' + config["server-config"]

if(config["rawdata-folder"][0]!='/'):
    config["rawdata-folder"] = config["project-folder"] + '/' + config["rawdata-folder"]

if(config["samplesheet-file"][0]!='/'):
    config["samplesheet-file"] = config["project-folder"] + '/' + config["samplesheet-file"]

if config["sampleinfo-file"] == "":
    pass
else:
    if(config["sampleinfo-file"][0]!='/'):
        config["sampleinfo-file"] = config["project-folder"] + '/' + config["sampleinfo-file"]

if(config["genome"][0]!='/'):
    config["genome"] = config["project-folder"] + '/' + config["genome"]

if(config["annot"][0]!='/'):
    config["annot"] = config["project-folder"] + '/' + config["annot"]

if config["contamination-folder"] == "":
    pass
else:
    if(config["contamination-folder"][0]!='/'):
        config["contamination-folder"] = config["project-folder"] + '/' + config["contamination-folder"]

if(config["tmp-folder"][0]!='/'):
    config["tmp-folder"] = config["project-folder"] + '/' + config["tmp-folder"]

##### load config and sample sheets #####

samplesheet = pd.read_table(config["samplesheet-file"]).set_index("rawsample", drop=False)
rawsamples=list(samplesheet.rawsample)
samples=list(set(list(samplesheet.sample_name)))
lane=list(samplesheet.lane)

# Get the basename fastq inputs
possible_ext = [".fastq", ".fq.gz", ".fastq.gz", ".fasta", ".fa", ".fa.gz", ".fasta.gz"]
ext = ".null"

reads1_tmp = list(samplesheet.read1)
reads1_trim = []
for r in reads1_tmp:
    for e in possible_ext:
        if r.endswith(e):
            addThis = r[:-len(e)]
            reads1_trim += [addThis]
            ext=e

reads2_tmp = list(samplesheet.read2)
reads2_trim = []
for r in reads2_tmp:
    for e in possible_ext:
        if r.endswith(e):
            addThis = r[:-len(e)]
            reads2_trim += [addThis] 
            ext=e

##### Extract the cluster resource requests from the server config #####
cluster=dict()
if os.path.exists(config["server-config"]):
    with open(config["server-config"]) as yml:
        cluster = yaml.load(yml, Loader=yaml.FullLoader)

workdir: config["project-folder"]

wildcard_constraints:
    rawsamples="|".join(rawsamples),
    samples="|".join(samples),
    reads1_trim="|".join(reads1_trim),
    reads2_trim="|".join(reads2_trim)
    
    
##### Config autofile #####

config["params"]["featurecounts"]["attributetype_merged"]="gene_id"
config["params"]["featurecounts"]["featuretype_merged"]="transcript"

##### Complete the input configuration

config["star-index"] = config["project-folder"]+"/Reference/STAR2.7.5a"
config["report-script"] = config["pipeline-folder"]+"/scripts/workflow-report.Rmd"

##### Function definitions #####

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


def get_fastq_for_concatenating_read1(wildcards):
    r1 = samplesheet.loc[samplesheet["sample_name"] == wildcards.samples]["read1"]
    path = config["rawdata-folder"] + "/"
    output = [path + x for x in r1]
    return output   

def get_fastq_for_concatenating_read2(wildcards):
    r1 = samplesheet.loc[samplesheet["sample_name"] == wildcards.samples]["read2"]
    path = config["rawdata-folder"] + "/"
    output = [path + x for x in r1]
    return output   


##### Singularity container #####
config["singularity"] = {}
config["singularity"]["fastp"] = "docker://fischuu/fastp:0.20.1"
config["singularity"]["rsem"] = "docker://fischuu/rsem:1.3.3"
config["singularity"]["star"] = "docker://fischuu/star:2.7.5a"
config["singularity"]["stringtie"] = "docker://fischuu/stringtie:2.1.2"
config["singularity"]["report"] = "docker://fischuu/r-gbs:3.6.3-0.1"

##### Print the welcome screen #####
print("#################################################################################")
print("##### Welcome to the RNAseq pipeline")
print("##### version: "+version)
print("#####")
print("##### Pipeline configuration")
print("##### --------------------------------")
print("##### project-folder       : "+config["project-folder"])
print("##### pipeline-folder      : "+config["pipeline-folder"])
print("##### pipeline-config      : "+config["pipeline-config"])
print("##### server-config        : "+config["server-config"])
print("##### report-script        : "+config["report-script"])
print("#####")
print("##### Singularity configuration")
print("##### --------------------------------")
print("##### fastp   : "+config["singularity"]["fastp"])
print("##### star      : "+config["singularity"]["star"])
print("##### rsem       : "+config["singularity"]["rsem"])
print("##### stringtie : "+config["singularity"]["stringtie"])
print("##### report : "+config["singularity"]["report"])
print("##### Runtime-configurations")
print("##### --------------------------------")
print("##### Genome                : "+ config["genome"])
print("##### Annotation            : "+ config["annot"])


##### run complete pipeline #####

rule all:
    input:
      # Step1-Preparations
        "%s/chrName.txt" % (config["star-index"]),
        "%s/Reference/RSEM.chrlist" % (config["project-folder"]),
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

rule preparations:
    input:
        "%s/chrName.txt" % (config["star-index"])

rule meta:
    input:
      # Step1-Preparations
        "%s/chrName.txt" % (config["star-index"]),
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
