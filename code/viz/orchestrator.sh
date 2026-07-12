#!/bin/bash
source /mnt/scratche/slow/ghlab/qiuchen/miniforge3/etc/profile.d/conda.sh && conda activate /mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq
export NUMBA_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 EMBED_NJOBS=8 MPLCONFIGDIR=/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl
cd /Users/meng01/qiuchen/project/hackathon/perturb-seq
L=/Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/logs
st() { echo "[$(date)] $1" >> $L/orchestrator.log; }
run_wave() {  # $1=tag $2=embed-cmd-args
  local tag=$1; shift
  st "$tag embed start"
  python "$@" >> $L/${tag}_embed.log 2>&1
  if [ $? -ne 0 ]; then st "$tag embed FAILED"; return 1; fi
  st "$tag plot start"
  python viz_outputs/code/run_plot.py $tag >> $L/${tag}_plot.log 2>&1
  if [ $? -eq 0 ]; then touch $L/${tag}.done; st "$tag DONE"; else st "$tag plot FAILED"; fi
}

# WAVE 1: D4 Rest+Stim8hr (rebuild cleanly with thread cap; uses run_embed via a wave1 config)
run_wave D4_Rest_Stim8hr viz_outputs/code/run_embed.py viz_outputs/code/wave1_config.json

# WAVE 2: D4 timecourse (needs D4_Stim48hr — already downloaded)
st "WAVE2 wait D4_Stim48hr"; for i in $(seq 1 2880); do [ -f $L/dl_D4_Stim48hr.done ] && break; sleep 30; done
run_wave D4_timecourse viz_outputs/code/run_embed.py viz_outputs/code/wave2_config.json

# WAVE 3: full atlas (needs all 10 downloads + GO flag)
st "WAVE3 wait all downloads + GO"
for i in $(seq 1 5760); do
  n=$(ls $L/dl_*.done 2>/dev/null | wc -l)
  [ "$n" -ge 10 ] && [ -f $L/wave3_GO ] && break; sleep 30
done
run_wave AllDonors_AllConditions viz_outputs/code/run_embed.py viz_outputs/code/wave3_config.json
st "ORCHESTRATOR COMPLETE"
