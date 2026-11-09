#!/bin/bash
# =============================================================================
# hpc/de_small.sbatch - init (then round 0 or, warm-started, round 1), the round-0 initialisation
# (then round 1), continue after the cap, reconfig of a running run (no round is submitted), or
# assemble of a DE-ITS run.
# --export RUN, WHAT (init | round0 | continue | reconfig | assemble), MODEL, DATA, WARM_POP,
# WARM_MAPS, PHENO_W, PHENO_S, FOLD (init only). After assemble: hpc/m2_chain.sh.
# =============================================================================
#SBATCH -J m2de_s
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 01:00:00
source Code/hpc/de_common.sh
if [ "$WHAT" = "init" ]; then
  rtask --what init --run "$RUN" --model "$MODEL" --data "${DATA:-real}" --warm-pop "${WARM_POP:-}" \
    --warm-maps "${WARM_MAPS:-}" --pheno-w "${PHENO_W:-}" --pheno-s "${PHENO_S:-}" --fold "${FOLD:-}" && bash Code/hpc/de_next.sh "$RUN"
elif [ "$WHAT" = "round0" ]; then
  rtask --what round0 --run "$RUN" && bash Code/hpc/de_next.sh "$RUN"
elif [ "$WHAT" = "continue" ]; then
  rtask --what continue --run "$RUN" && bash Code/hpc/de_next.sh "$RUN"
elif [ "$WHAT" = "reconfig" ]; then
  rtask --what reconfig --run "$RUN"
else
  rtask --what assemble --run "$RUN" && bash Code/hpc/m2_chain.sh "$(dirname "$RUN")"
fi
