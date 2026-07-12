#!/bin/bash
# args: shard cond out_h5 out_sum seed coalesce log done
source /mnt/scratche/slow/ghlab/qiuchen/miniforge3/etc/profile.d/conda.sh
conda activate /mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq
export PYTHONUNBUFFERED=1
python /Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_scripts/phaseE_gather_targeted.py "$1" "$2" "$3" "$4" "$5" "$6" > "$7" 2>&1
echo "EXIT=$?" > "$8"
