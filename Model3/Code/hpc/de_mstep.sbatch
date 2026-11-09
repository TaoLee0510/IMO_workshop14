#!/bin/bash
# =============================================================================
# hpc/de_mstep.sbatch - DE-ITS M step (closed forms + structural DEoptim); afterwards
# hpc/de_next.sh reads the status and submits the next round or the final E step.
# --export RUN, ITER.
# =============================================================================
#SBATCH -J m2de_m
#SBATCH -c 60
#SBATCH --mem=96G
#SBATCH -t 11:00:00
source Code/hpc/de_common.sh
rtask --what mstep --run "$RUN" --iter "$ITER" --cores "$SLURM_CPUS_PER_TASK"
bash Code/hpc/de_next.sh "$RUN"
