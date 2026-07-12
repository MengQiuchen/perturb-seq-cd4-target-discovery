#!/bin/bash
cd /Users/meng01/qiuchen/project/hackathon/perturb-seq
export MPLCONFIGDIR=/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl
/mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq/bin/python phaseE_scripts/diag_read.py > /Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/diag_read.log 2>&1
echo "EXIT=$?" > /Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/diag_read.done
