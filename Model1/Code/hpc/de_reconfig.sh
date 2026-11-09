#!/bin/bash
# =============================================================================
# hpc/de_reconfig.sh RUN... - move running DE-ITS runs to the current outer-loop settings
# (Code/config.R) from their next M step on: one small job per run rewrites meta.rds (the
# earlier settings and the round are kept in meta$reconfigured); nothing is submitted.
# Run from the task folder.
# =============================================================================
source Code/hpc/de_common.sh
for RUN in "$@"; do
  j=$(sub -J "r_$(basename "$RUN")" -o "$RUN/logs/reconfig_%j.out" \
      --export=ALL,RUN="$RUN",WHAT=reconfig Code/hpc/de_small.sbatch)
  echo "$(date '+%F %T') reconfig job $j" >> "$RUN/submissions.log"
  echo "$RUN reconfig $j"
done
