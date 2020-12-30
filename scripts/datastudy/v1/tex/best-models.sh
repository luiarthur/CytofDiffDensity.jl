#!/bin/bash

function grep_dic {
  info_dir=$1
  cat $info_dir/info.txt | grep -oP '(?<=curves:\s)\d+\.\d+(?=\s\(exponentiate=t)'
}


export best=(\
  "results/K=3_marker=CD3z_skewtmix=true/"
  "results/K=4_marker=CD3z_skewtmix=false/"

  "results/K=5_marker=CD56_skewtmix=true/"
  "results/K=5_marker=CD56_skewtmix=false/"

  "results/K=3_marker=CD57_skewtmix=true/"
  "results/K=2_marker=CD57_skewtmix=false/"

  "results/K=3_marker=EOMES_skewtmix=true/"
  "results/K=4_marker=EOMES_skewtmix=false/"

  "results/K=5_marker=Granzyme_A_skewtmix=true/"
  "results/K=5_marker=Granzyme_A_skewtmix=false/"

  "results/K=4_marker=LAG3_skewtmix=true/"
  "results/K=4_marker=LAG3_skewtmix=false/"

  "results/K=3_marker=Perforin_skewtmix=true/"
  "results/K=4_marker=Perforin_skewtmix=false/"

  "results/K=3_marker=Siglec7_skewtmix=true/"
  "results/K=6_marker=Siglec7_skewtmix=false/"
)

for b in ${best[@]}
do
  echo $b `grep_dic $b`
done
