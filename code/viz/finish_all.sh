#!/bin/bash
# Consolidated autonomous finisher: completes Waves 2 & 3 and packages everything.
# Idempotent — skips any wave whose .done marker already exists.
# Waves 1-2 run automatically; Wave 3 (untested cross-donor Harmony atlas) is GATED on an
# explicit `touch viz_outputs/logs/wave3_GO` flag — see the WAVE 3 block below.
source /mnt/scratche/slow/ghlab/qiuchen/miniforge3/etc/profile.d/conda.sh
conda activate /mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq
export NUMBA_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 EMBED_NJOBS=8
export MPLCONFIGDIR=/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl
PROJ=/Users/meng01/qiuchen/project/hackathon/perturb-seq
cd $PROJ
L=$PROJ/viz_outputs/logs
st(){ echo "[$(date)] $1" >> $L/finish_all.log; }

run_wave(){  # $1=tag $2=config
  local tag=$1 cfg=$2
  if [ -f $L/${tag}.done ]; then st "$tag already done, skip"; return 0; fi
  st "$tag embed start"
  python viz_outputs/code/run_embed.py $cfg >> $L/${tag}_embed.log 2>&1
  if [ $? -ne 0 ]; then st "$tag embed FAILED"; return 1; fi
  st "$tag plot start"
  python viz_outputs/code/run_plot.py $tag >> $L/${tag}_plot.log 2>&1
  if [ $? -eq 0 ]; then touch $L/${tag}.done; st "$tag DONE"; else st "$tag plot FAILED"; return 1; fi
}

package(){  # rebuild manifest + figures tarball from whatever waves are done so far
  st "packaging: manifest + tar"
  python viz_outputs/code/make_manifest.py >> $L/finish_all.log 2>&1
  ( cd $PROJ/viz_outputs && tar czf figures_bundle.tar.gz figures/ FIGURES_README.md 2>> $L/finish_all.log )
}

st "FINISH_ALL START"

# WAVE 1 (should already be done; rebuild if marker missing but checkpoint exists it will reuse cache fast)
run_wave D4_Rest_Stim8hr viz_outputs/code/wave1_config.json

# WAVE 2: D4 timecourse — needs D4_Stim48hr
st "WAVE2 wait D4_Stim48hr download"
for i in $(seq 1 2880); do [ -f $L/dl_D4_Stim48hr.done ] && break; sleep 30; done
run_wave D4_timecourse viz_outputs/code/wave2_config.json

# Package Waves 1+2 now so their figures are delivered regardless of the Wave 3 gate.
package
touch $L/waves12_packaged.done
st "WAVES 1+2 packaged"

# WAVE 3: full atlas w/ Harmony — needs all 10 downloads AND explicit GO flag.
# GATE RETAINED: Wave 3 is the untested path (Harmony integration over donor across all
# 12 shards, never yet exercised). Wave 1 validation does NOT extend to it. Requires a
# deliberate `touch $L/wave3_GO` after reviewing the D4 figures. Do not remove this gate.
st "WAVE3 wait for all 10 downloads + wave3_GO flag"
for i in $(seq 1 5760); do
  n=$(ls $L/dl_*.done 2>/dev/null | wc -l)
  [ "$n" -ge 10 ] && [ -f $L/wave3_GO ] && break; sleep 30
done
st "WAVE3 gate open (downloads=$(ls $L/dl_*.done 2>/dev/null | wc -l), GO=$([ -f $L/wave3_GO ] && echo yes || echo no))"
run_wave AllDonors_AllConditions viz_outputs/code/wave3_config.json

# Final package (now includes Wave 3 figures if the gate opened and it ran)
package
touch $L/finish_all.done
st "FINISH_ALL COMPLETE"
