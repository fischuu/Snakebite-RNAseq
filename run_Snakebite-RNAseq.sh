module load snakemake

pipelineFolder="/users/fischerd/git/Snakebite-RNAseq"
projectFolder="/scratch/project_2001746/RNAseq_Example"

export APPTAINER_TMPDIR="/scratch/project_2001746/tmp"
export APPTAINER_CACHEDIR="/scratch/project_2001746/tmp"

mkdir -p $APPTAINER_TMPDIR

# Create the rulegraph
snakemake -s $pipelineFolder/Snakebite-RNAseq.smk \
          --configfile $projectFolder/Snakebite-RNAseq_config.yaml \
          --rulegraph | dot -T png > $projectFolder/workflow.png

snakemake -s $pipelineFolder/Snakebite-RNAseq.smk \
          -j 300 \
          --use-singularity \
          --singularity-args "-B /scratch:/scratch" \
          --configfile $projectFolder/Snakebite-RNAseq_config.yaml \
          --latency-wait 60 \
          --cluster-config $projectFolder/Snakebite-RNAseq_server_config.yaml \
          --cluster "sbatch -t {cluster.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} -p {cluster.partition} -D {cluster.working-directory}" \
          --scheduler greedy \
	--cluster-cancel scancel \
          $@
