#!/bin/bash
#SBATCH --job-name=pt_umap_full
#SBATCH --partition=cuda
#SBATCH --qos=cuda_limit
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/slurm_pt_full_%j.out
#SBATCH --error=/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/slurm_pt_full_%j.out

set -e
PROJ=/Users/meng01/qiuchen/project/hackathon/perturb-seq
PY=/mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq/bin/python
OUT=$PROJ/phaseE_outputs/umap_within

echo "=== JOB $SLURM_JOB_ID on $(hostname) at $(date) ==="
echo "=== node: $SLURMD_NODENAME  cpus: $SLURM_CPUS_PER_TASK ==="
nvidia-smi -L 2>/dev/null || echo "(no nvidia-smi / GPU not used — CPU scanpy)"

# I/O sanity check on this node's mount to the raw shards (pool 'slow')
echo "=== I/O probe: 256MB fresh read from D4_Rest shard ==="
dd if=$PROJ/perturb-seq_data/cell_level/D4_Rest.assigned_guide.h5ad of=/dev/null \
   bs=8M count=32 iflag=skip_bytes skip=40000000000 2>&1 | tail -1 || true

cd $PROJ
export MPLCONFIGDIR=/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl
export EMBED_NJOBS=16
export OMP_NUM_THREADS=16

echo "=== starting per_target_umap.py (full uncapped gather) at $(date) ==="
$PY phaseE_scripts/per_target_umap.py
echo "=== EXIT=$? at $(date) ==="
