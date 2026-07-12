#!/bin/bash
source /mnt/scratche/slow/ghlab/qiuchen/miniforge3/etc/profile.d/conda.sh
conda activate /mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq
export PYTHONUNBUFFERED=1
python /Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_scripts/phaseE_analyze.py "$1" "$2" "$3" "$4" "$5" > "$6" 2>&1
echo "EXIT=$?" > "$7"
