#!/bin/bash
source /mnt/scratche/slow/ghlab/qiuchen/miniforge3/etc/profile.d/conda.sh && conda activate /mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq
export MPLCONFIGDIR=/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl
cd /Users/meng01/qiuchen/project/hackathon/perturb-seq
# wait for embed to finish (marker = "DONE wave1" in log)
for i in $(seq 1 240); do
  if grep -q "DONE wave1" /Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/logs/wave1_embed.log 2>/dev/null && [ -f /Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/checkpoints/D4_Rest_Stim8hr.embedded.h5ad ]; then break; fi
  sleep 15
done
echo "PLOT_START $(date)" >> /Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/logs/wave1_plot.log
python viz_outputs/code/wave1_plot.py >> /Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/logs/wave1_plot.log 2>&1
rc=$?
echo "PLOT_DONE rc=$rc $(date)" >> /Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/logs/wave1_plot.log
if [ $rc -eq 0 ]; then touch /Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/logs/wave1_plot.done; fi
