#!/bin/bash

# Run this script with: sbatch submit-job.sh

#SBATCH -p 128x24      # Partition name
#SBATCH -J data-v4 # Job name
#SBATCH --mail-user=alui2@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o out/slurm-job.out    # Name of stdout output file
#SBATCH -N 3         # Total number of nodes requested (128x24/Instructional only)
#SBATCH -t 120:00:00  # Run Time (hh:mm:ss)
#SBATCH --mem=42G # Memory to be allocated PER NODE

date 
echo "nproc: `nproc`"
echo "Starting jobs!"
julia run.jl &
echo "Done submitting jobs."

echo "Job submission time:"
date

echo "Jobs are now running. A message will be printed and emailed when jobs are done."
wait

echo "Jobs are completed."
echo "Job completion time:"
date
