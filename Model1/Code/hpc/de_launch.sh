#!/bin/bash
# =============================================================================
# hpc/de_launch.sh RUN MODEL [DATA] - start a Model1 DE-ITS fit: submits the init job, which
# writes the run directory and submits round 0 (or round 1 when warm-started); every M step then submits the
# next round itself (hpc/de_next.sh) until convergence, the final E step and assemble.
# Run from the task folder that holds Code/, Data/derived/ and Results/.
# DATA: real (default; requires the Model1 pre-registration) or a synthetic cohort .rds.
# WARM_POP / WARM_MAPS (environment, optional): start from another run's population file and
# per-patient E-step directory (R/fit_de.R de_init). PHENO_W / PHENO_S (optional): weight and smoothing
# scale of the phenotype term (S2: w 0; sensitivity analyses). FOLD (optional): V1 training set without fold FOLD.
# =============================================================================
RUN="$1"; MODEL="$2"; DATA="${3:-real}"
source Code/hpc/de_common.sh
mkdir -p "$RUN/logs"
j=$(sub -J "i_$(basename "$RUN")" -o "$RUN/logs/init_%j.out" \
    --export=ALL,RUN="$RUN",WHAT=init,MODEL="$MODEL",DATA="$DATA",WARM_POP="${WARM_POP:-}",WARM_MAPS="${WARM_MAPS:-}",PHENO_W="${PHENO_W:-}",PHENO_S="${PHENO_S:-}",FOLD="${FOLD:-}" \
    Code/hpc/de_small.sbatch)
echo "$(date '+%F %T') init job $j model $MODEL data $DATA warm ${WARM_POP:-none} pheno_w ${PHENO_W:-default} pheno_s ${PHENO_S:-default} fold ${FOLD:-none}" >> "$RUN/submissions.log"
echo "$j"
