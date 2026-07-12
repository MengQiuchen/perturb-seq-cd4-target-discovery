#!/bin/bash
#SBATCH --job-name=xdonor_within
#SBATCH --partition=cuda
#SBATCH --qos=cuda_limit
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/slurm_xdonor_%j.out
#SBATCH --error=/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/slurm_xdonor_%j.out

set -e
PROJ=/Users/meng01/qiuchen/project/hackathon/perturb-seq
PY=/mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq/bin/python

echo "=== JOB $SLURM_JOB_ID on $(hostname) at $(date) ==="
echo "=== node: $SLURMD_NODENAME cpus: $SLURM_CPUS_PER_TASK ==="

cd $PROJ
export MPLCONFIGDIR=/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl
export EMBED_NJOBS=16
export OMP_NUM_THREADS=16

echo "=== starting cross_donor_within.py at $(date) ==="
$PY phaseE_scripts/cross_donor_within.py
echo "=== EXIT=$? at $(date) ==="
