# Pipeline-RNA-seq
Standard RNA-seq snakemake pipeline.

## Installation
The test data is not part of this repository and needs to be downloaded separately here:
 
```
wget ...
```

To clone into the pipeline type

```
git clone https://github.com/fischuu/Pipeline-RNA-seq.git
```
## Usage
The pipeline is configured in the file `...` with the following required parameters.

In addition, the follwoing must be provided

## Preparations:

### Create rawsamples files

### Create samples files

### Create sampleInfo.txt file
cut -d'-' -f1 samples > animal
cut -d'-' -f2 samples > pathogen
cut -d'-' -f3 samples > time
cut -d'-' -f4 samples > replicate
cut -d'-' -f5 samples > sequencer


paste -d'\t' samples animal pathogen time replicate sequencer > sampleList.tmp

rm animal
rm time
rm sequencer
rm pathogen
rm replicate

printf "sample\tanimal\tpathogen\ttime\treplicate\tsequencer\n" > header
cat header sampleList.tmp > sampleInfo.txt
rm sampleList.tmp
rm header

# Current status
Under development and possible not yet stable for production
