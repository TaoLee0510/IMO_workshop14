#!/bin/bash
# =============================================================================
# hpc/de_small.sbatch - init (then the first round), continue after the cap, reconfig of a running
# run (no round is submitted), or assemble of a DE-ITS run.
# --export RUN, WHAT (init | continue | reconfig | assemble), MODEL, DATA, WARM_POP, WARM_MAPS (init only).
# =============================================================================
#SBATCH -J a4de_s
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 01:00:00
source Code/hpc/de_common.sh
if [ "$WHAT" = "init" ]; then
  rtask --what init --run "$RUN" --model "$MODEL" --data "${DATA:-real}" --warm-pop "${WARM_POP:-}" \
    --warm-maps "${WARM_MAPS:-}" && bash Code/hpc/de_next.sh "$RUN"
elif [ "$WHAT" = "continue" ]; then
  rtask --what continue --run "$RUN" && bash Code/hpc/de_next.sh "$RUN"
elif [ "$WHAT" = "reconfig" ]; then
  rtask --what reconfig --run "$RUN"
else
  rtask --what assemble --run "$RUN"
fi
