#!/bin/bash

# Run this script with: sbatch submit-job.sh

#SBATCH -p 128x24      # Partition name
#SBATCH -J post-v4 # Job name
#SBATCH --mail-user=alui2@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o out/slurm-post.out    # Name of stdout output file
#SBATCH -N 1  # number of nodes
#SBATCH -t 5:00:00  # Run Time (hh:mm:ss)
#SBATCH --mem=42G # Memory to be allocated PER NODE

date 
echo "nproc: `nproc`"
echo "SLURM_JOB_NUM_NODES: $SLURM_JOB_NUM_NODES"
echo "SLURM_CPUS_ON_NODE: $SLURM_CPUS_ON_NODE"

echo "Job submission time:"
date

echo "Concatenating results ..."
julia concat-results.jl
echo "Done"

echo "Jobs are completed."
echo "Job completion time:"
date
