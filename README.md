# Pipeline-RNA-seq
Standard RNA-seq pipeline.

## Installation
The test data is not part of this repository and needs to be downloaded separately here:
 
```
wget ...
```

To clone into the pipeline type

```
git clone git@github.com:fischuu/Pipeline-RNA-seq.git
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
cut -d'-' -f4 samples > sequencer


paste -d'\t' samples animal pathogen time sequencer > tmp

rm animal
rm time
rm sequencer
rm pathogen

printf "sample\tanimal\tpathogen\ttime\tsequencer\n" > header
cat header tmp > sampleInfo.txt
rm tmp
rm header
