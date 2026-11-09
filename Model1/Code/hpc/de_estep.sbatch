#!/bin/bash
# =============================================================================
# hpc/de_estep.sbatch - DE-ITS E step (or final E step with the importance-sampling
# marginal likelihood), one array task per patient. Submitted by hpc/de_next.sh with
# --export RUN, ITER, WHAT (estep | efinal) and --array=1-n.
# =============================================================================
#SBATCH -J m1de_e
#SBATCH -c 22
#SBATCH --mem=24G
#SBATCH -t 06:00:00
source Code/hpc/de_common.sh
rtask --what "$WHAT" --run "$RUN" --iter "$ITER" --patient "$SLURM_ARRAY_TASK_ID" --cores "$SLURM_CPUS_PER_TASK"
