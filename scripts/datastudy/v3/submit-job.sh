#!/bin/bash

# Run this script with: sbatch submit-job.sh

#SBATCH -p 128x24      # Partition name
#SBATCH -J data-v4 # Job name
#SBATCH --mail-user=alui2@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o out/slurm-job.out    # Name of stdout output file
#SBATCH -n 128        # number of tasks
#SBATCH -c 1         # Number of cores per task
#SBATCH -t 120:00:00  # Run Time (hh:mm:ss)
#SBATCH --mem=42G # Memory to be allocated PER NODE

date 
echo "nproc: `nproc`"
echo "SLURM_JOB_NUM_NODES: $SLURM_JOB_NUM_NODES"
echo "SLURM_CPUS_ON_NODE: $SLURM_CPUS_ON_NODE"

# TODO: Rewrite this. One task can only use up to the max number of cpus on the node.
# 1 task cannot be split into multiple nodes.

# echo "Starting jobs!"
# julia run.jl &
# echo "Done submitting jobs."
# 
# echo "Job submission time:"
# date
# 
# echo "Jobs are now running. A message will be printed and emailed when jobs are done."
# wait
# 
# echo "Jobs are completed."
# echo "Job completion time:"
# date
