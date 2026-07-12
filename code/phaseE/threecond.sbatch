#!/bin/bash
#SBATCH --job-name=threecond
#SBATCH --partition=cuda
#SBATCH --qos=cuda_limit
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/slurm_threecond_%j.out
#SBATCH --error=/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/slurm_threecond_%j.out

PROJ=/Users/meng01/qiuchen/project/hackathon/perturb-seq
PY=/mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq/bin/python
echo "=== JOB $SLURM_JOB_ID on $SLURMD_NODENAME at $(date) ==="
cd $PROJ
export MPLCONFIGDIR=/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl
export EMBED_NJOBS=16
export OMP_NUM_THREADS=16
$PY phaseE_scripts/three_cond_within.py
echo "=== EXIT=$? at $(date) ==="
