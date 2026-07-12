#!/bin/bash
L=/Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/logs
echo "### $(date) ###"
echo "-- procs --"; pgrep -af 'run_embed|wave1_run|run_plot|wave1_plot|aws s3 cp|_watch' | grep -v pgrep | sed 's/ .*bin\///' | cut -c1-90
echo "-- wave1 embed --"; tail -2 $L/wave1_embed.log 2>/dev/null
echo "-- wave1 plot --"; tail -3 $L/wave1_plot.log 2>/dev/null
echo "-- downloads done --"; ls $L/dl_*.done 2>/dev/null | sed 's#.*/dl_##;s/.done//' | tr '\n' ' '; echo
echo "-- download progress --"; grep -E 'DONE|COMPLETE' $L/dl_master.log 2>/dev/null | tail -4
echo "-- figures --"; ls /Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/figures/*.png 2>/dev/null | wc -l
echo "-- shardcache --"; ls /Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/checkpoints/shardcache/*.sub.h5ad 2>/dev/null | wc -l
