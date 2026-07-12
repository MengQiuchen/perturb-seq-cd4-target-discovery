#!/bin/bash
cd /Users/meng01/qiuchen/project/hackathon/perturb-seq
export MPLCONFIGDIR=/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl
export EMBED_NJOBS=8
/mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq/bin/python phaseE_scripts/per_target_subset.py > /Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/pt_subset.log 2>&1
echo "EXIT=$?" > /Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/pt_subset.done
