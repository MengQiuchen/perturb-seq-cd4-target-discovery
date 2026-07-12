#!/bin/bash
set -uo pipefail
source /mnt/scratche/slow/ghlab/qiuchen/miniforge3/etc/profile.d/conda.sh
conda activate /mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq
export NUMBA_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 EMBED_NJOBS=8
export MPLCONFIGDIR=/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl
PROJ=/Users/meng01/qiuchen/project/hackathon/perturb-seq
cd $PROJ
L=$PROJ/viz_outputs/logs
CODE=$PROJ/viz_outputs/code
CK=$PROJ/viz_outputs/checkpoints
st(){ echo "[$(date '+%F %T')] $1" >> $L/backfill.log; }

run_donor(){  # $1=tag $2=config
  local tag=$1 cfg=$2
  if [ -f $L/${tag}.done ]; then st "$tag already done, skip"; return 0; fi
  st "$tag EMBED start"
  python $CODE/run_embed.py $cfg >> $L/${tag}_embed.log 2>&1
  if [ $? -ne 0 ]; then st "$tag EMBED FAILED (see ${tag}_embed.log)"; return 1; fi
  st "$tag PLOT start"
  python $CODE/run_plot.py $tag >> $L/${tag}_plot.log 2>&1 || { st "$tag PLOT FAILED"; return 1; }
  local ckpt=$CK/${tag}.embedded.h5ad
  st "$tag BYCOND start"
  python $CODE/per_target_bycondition.py $ckpt $tag >> $L/${tag}_plot.log 2>&1 || st "$tag per_target_bycondition WARN"
  python $CODE/per_module_bycondition.py $ckpt $tag >> $L/${tag}_plot.log 2>&1 || st "$tag per_module_bycondition WARN"
  touch $L/${tag}.done; st "$tag DONE ($(ls $L/../figures/*${tag}*.png 2>/dev/null | wc -l) figures)"
}

st "BACKFILL START (D1/D2/D3 per-donor timecourse)"
run_donor D1_timecourse $CODE/d1_timecourse_config.json
run_donor D2_timecourse $CODE/d2_timecourse_config.json
run_donor D3_timecourse $CODE/d3_timecourse_config.json
st "rebuild manifest + tarball"
python $CODE/make_manifest.py >> $L/backfill.log 2>&1
( cd $PROJ/viz_outputs && tar czf figures_bundle.tar.gz figures/ FIGURES_README.md 2>> $L/backfill.log )
touch $L/backfill_all.done
st "BACKFILL COMPLETE"
