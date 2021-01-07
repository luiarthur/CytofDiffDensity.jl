#!/bin/bash

jid1=$(sbatch slurm-spawn-jobs.sh | grep -oP '\d+'); echo $jid1
jid2=$(sbatch --dependency=afterok:$jid1 slurm-concat-results.sh | grep -oP '\d+'); echo $jid2
