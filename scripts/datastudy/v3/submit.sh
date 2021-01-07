#!/bin/bash

jid1=$(sbatch spawn-jobs.sh | grep -oP '\d+'); echo $jid1
jid2=$(sbatch --dependency=afterok:$jid1 concat-results.sh | grep -oP '\d+'); echo $jid2
