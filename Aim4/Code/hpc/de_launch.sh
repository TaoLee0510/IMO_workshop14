#!/bin/bash
# =============================================================================
# hpc/de_launch.sh RUN MODEL [DATA] - start a DE-ITS fit (v1.3): submits the init job,
# which writes the run directory and submits round 1; every M step then submits the
# next round itself (hpc/de_next.sh) until convergence, the final E step and assemble.
# Run from the task folder that holds Code/, Data/derived_v13/ and Results/.
# DATA: real (default; requires the v1.3 pre-registration) or a synthetic cohort .rds.
# WARM_POP / WARM_MAPS (environment, optional): start from another run's population file and
# per-patient E-step directory (R/fit_de.R de_init).
# =============================================================================
RUN="$1"; MODEL="$2"; DATA="${3:-real}"
source Code/hpc/de_common.sh
mkdir -p "$RUN/logs"
j=$(sub -J "i_$(basename "$RUN")" -o "$RUN/logs/init_%j.out" \
    --export=ALL,RUN="$RUN",WHAT=init,MODEL="$MODEL",DATA="$DATA",WARM_POP="${WARM_POP:-}",WARM_MAPS="${WARM_MAPS:-}" \
    Code/hpc/de_small.sbatch)
echo "$(date '+%F %T') init job $j model $MODEL data $DATA warm ${WARM_POP:-none}" >> "$RUN/submissions.log"
echo "$j"
