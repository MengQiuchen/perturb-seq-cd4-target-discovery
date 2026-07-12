#!/bin/bash
cd /Users/meng01/qiuchen/project/hackathon/perturb-seq
export MPLCONFIGDIR=/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl
/mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq/bin/python phaseE_scripts/inspect_cache.py > /Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/inspect_cache.log 2>&1
echo "EXIT=$?" > /Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/inspect_cache.done
