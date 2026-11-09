#!/bin/bash
# =============================================================================
# hpc/de_continue.sh RUN... - continue DE-ITS runs that stopped at the outer safety cap under
# the current outer-loop settings (Code/config.R): one small job per run keeps the cap's final
# E step and fit (final_capNNN, fit_capNNN.rds), sets the status to the round after the cap
# and submits that round (hpc/de_next.sh). Run from the task folder.
# =============================================================================
source Code/hpc/de_common.sh
for RUN in "$@"; do
  j=$(sub -J "i_$(basename "$RUN")" -o "$RUN/logs/continue_%j.out" \
      --export=ALL,RUN="$RUN",WHAT=continue Code/hpc/de_small.sbatch)
  echo "$(date '+%F %T') continue job $j" >> "$RUN/submissions.log"
  echo "$RUN continue $j"
done
