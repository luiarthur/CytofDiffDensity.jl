#!/bin/bash

# Run this script with: sbatch submit-job.sh

#SBATCH -p 128x24      # Partition name
#SBATCH -J data-v3 # Job name
#SBATCH --mail-user=alui2@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o out/slurm-job.out    # Name of stdout output file
#SBATCH -N 3  # number of nodes
#SBATCH --ntasks-per-node 24  # maximum number of tasks to be run on a node.
#SBATCH -t 120:00:00  # Run Time (hh:mm:ss)
#SBATCH --mem=42G # Memory to be allocated PER NODE

date 
echo "nproc: `nproc`"
echo "SLURM_JOB_NUM_NODES: $SLURM_JOB_NUM_NODES"
echo "SLURM_CPUS_ON_NODE: $SLURM_CPUS_ON_NODE"

# NOTE: One task can only use up to the max number of cpus on the node.
# One task cannot be split into multiple nodes.

markers="CD3z EOMES Perforin Granzyme_A Siglec7 LAG3 CD56 CD57"
Ks=`seq 2 9`
skewmixt="0 1"

echo "Compiling job ..."
srun -n1 -c1 julia imports.jl

echo "Starting jobs!"
for marker in $markers
do
  for K in $Ks
  do
    for stm in $skewmixt
    do
      logname="marker=${marker}_K=${K}_stm=${stm}"
      srun -n1 -c1 julia run.jl $marker $K $stm &> out/${logname}-log.txt &
    done
  done
done
echo "Done submitting jobs."

echo "Job submission time:"
date

echo "Jobs are now running. A message will be printed and emailed when jobs are done."
wait

echo "Jobs are completed."
echo "Job completion time:"
date
