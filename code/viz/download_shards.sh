#!/bin/bash
set -u
D=/Users/meng01/qiuchen/project/hackathon/perturb-seq/perturb-seq_data/cell_level
LOG=/Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/logs
BUCKET=s3://genome-scale-tcell-perturb-seq/marson2025_data
source /mnt/scratche/slow/ghlab/qiuchen/miniforge3/etc/profile.d/conda.sh
conda activate /mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq
# priority order: D4_Stim48hr first (Wave 2), then D1/D2/D3 all conditions (Wave 3)
FILES=(D4_Stim48hr D1_Rest D1_Stim8hr D1_Stim48hr D2_Rest D2_Stim8hr D2_Stim48hr D3_Rest D3_Stim8hr D3_Stim48hr)
MAXJOBS=3
download_one() {
  local f=$1
  local dest=$D/${f}.assigned_guide.h5ad
  if [ -f $LOG/dl_${f}.done ]; then echo "SKIP $f (already done)" >> $LOG/dl_master.log; return; fi
  local start=$(date +%s)
  echo "START $f $(date)" >> $LOG/dl_master.log
  aws s3 cp --no-sign-request $BUCKET/${f}.assigned_guide.h5ad $dest > $LOG/dl_${f}.log 2>&1
  local rc=$?
  local end=$(date +%s)
  if [ $rc -eq 0 ]; then touch $LOG/dl_${f}.done; fi
  echo "DONE $f rc=$rc secs=$((end-start)) $(date)" >> $LOG/dl_master.log
}
echo "===== DOWNLOAD RUN START $(date) =====" >> $LOG/dl_master.log
for f in "${FILES[@]}"; do
  while [ $(jobs -r | wc -l) -ge $MAXJOBS ]; do sleep 10; done
  download_one "$f" &
done
wait
echo "===== ALL_DOWNLOADS_COMPLETE $(date) =====" >> $LOG/dl_master.log
