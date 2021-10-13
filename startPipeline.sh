# Run the RNA-Seq snakemake pipeline on the puhti cluster
#
# Before you start, take care that the following files and folder are available:

# STAR Index (folder)
################################################################################
mkdir -p /scratch/project_2002561/MastitisChallenge/RNASeq-Analysis/Reference/STAR2.7.5a

# SLURM logs (folder)
################################################################################
mkdir -p /scratch/project_2002561/MastitisChallenge/logs/

# rawsamples (file)
###########################
find /scratch/project_2002561/MastitisChallenge/RNASeq-Analysis/FASTQ/RAW/ -name '*R1_001.fastq.gz' | xargs -n1 basename | sed 's/_R1_001.fastq.gz//g' >  /scratch/project_2002561/MastitisChallenge/RNASeq-Analysis/rawsamples

# samples (file)
###########################
ls /scratch/project_2002561/MastitisChallenge/RNASeq-Analysis/FASTQ/RAW/ | cut -d '_' -f1 | sort | uniq > /scratch/project_2002561/MastitisChallenge/RNASeq-Analysis/samples

module load bioconda/3
source activate LukeNGS

# Create the rulegraph
#snakemake -s /scratch/project_2002561/pipeline/Snakefile-Pipeline-RNA-seq.smk \
#          --configfile /scratch/project_2002561/MastitisChallenge/RNASeq-Analysis/Pipeline-RNA-seq_config.yaml \
#          --rulegraph | dot -T png > /scratch/project_2002561/MastitisChallenge/RNASeq-Analysis/workflow.png

snakemake -s /scratch/project_2002561/pipeline/Snakefile-Pipeline-RNA-seq.smk \
          -j 300 \
          --use-singularity \
          --singularity-args "-B /scratch:/scratch" \
          --configfile /scratch/project_2002561/MastitisChallenge/RNASeq-Analysis/Pipeline-RNA-seq_config.yaml \
          --latency-wait 60 \
          --cluster-config /scratch/project_2002561/pipeline/Pipeline-RNA-seq_server_config.yaml \
          --cluster "sbatch -t {cluster.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} -p {cluster.partition} -D {cluster.working-directory}" \
          $1
