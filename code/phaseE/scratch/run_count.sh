#!/bin/bash
cd /Users/meng01/qiuchen/project/hackathon/perturb-seq
export MPLCONFIGDIR=/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl
mkdir -p /Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within
/mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq/bin/python phaseE_scripts/count_obs.py > /Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/count_obs.log 2>&1
echo "EXIT=$?" > /Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/count_obs.done
