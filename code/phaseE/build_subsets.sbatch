#!/bin/bash
#SBATCH --job-name=build_subsets
#SBATCH --partition=cuda
#SBATCH --qos=cuda_limit
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --time=06:00:00
#SBATCH --output=/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/slurm_build_%j.out
#SBATCH --error=/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/slurm_build_%j.out

PROJ=/Users/meng01/qiuchen/project/hackathon/perturb-seq
PY=/mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq/bin/python
CELL=$PROJ/perturb-seq_data/cell_level
CKD=$PROJ/phaseE_outputs/checkpoints
GATHER=$PROJ/phaseE_scripts/phaseE_gather.py

echo "=== JOB $SLURM_JOB_ID on $SLURMD_NODENAME at $(date) ==="
cd $PROJ
export MPLCONFIGDIR=/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl

# 7 missing donor x condition subsets (seed=1 to match existing Stim8hr subsets)
# format: DONOR COND
COMBOS="D1 Rest
D2 Rest
D3 Rest
D4 Rest
D1 Stim48hr
D2 Stim48hr
D3 Stim48hr"

echo "$COMBOS" | while read D COND; do
  [ -z "$D" ] && continue
  OUT=$CKD/${D}_${COND}.subset.h5ad
  SUM=$CKD/${D}_${COND}.subset.summary.json
  SHARD=$CELL/${D}_${COND}.assigned_guide.h5ad
  if [ -f "$OUT" ]; then
    echo "=== SKIP ${D}_${COND} (exists) ==="
    continue
  fi
  echo "=== BUILD ${D}_${COND} from $SHARD at $(date) ==="
  $PY $GATHER "$SHARD" "$COND" "$OUT" "$SUM" 1
  echo "=== DONE ${D}_${COND} rc=$? at $(date) ==="
done

echo "=== ALL BUILDS DONE at $(date) ==="
