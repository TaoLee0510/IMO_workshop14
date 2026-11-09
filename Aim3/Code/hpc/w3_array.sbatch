#!/bin/bash
# =============================================================================
# hpc/w3_array.sbatch
# Purpose : W3 on the cluster, one array task per (scenario x replicate x n):
#           4 scenarios x 20 replicates at n = 35 (tasks 1-80) and
#           4 scenarios x 2 replicates at n = 350 (tasks 81-88).
# Usage   : from the Aim3 folder on the cluster (holds Code/ and
#           Results/real/objects/W0.rds; Rlib/ with Rcpp one level up):
#             mkdir -p logs && sbatch Code/hpc/w3_array.sbatch
#           Then copy Results/real/objects/W3_parts/ back and run
#             Rscript Code/run_all.R --steps W3
# =============================================================================
#SBATCH -J aim3_w3
#SBATCH -p red
#SBATCH --qos=xxlarge
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH -t 10:00:00
#SBATCH -a 1-88
#SBATCH -o logs/w3_%A_%a.out

IMG=${IMG:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif}
RLIB=${RLIB:-$PWD/../Rlib}
export SINGULARITYENV_R_LIBS="$RLIB"
export APPTAINERENV_R_LIBS="$RLIB"
singularity exec --cleanenv -B /share,/tmp "$IMG" \
  Rscript Code/hpc/w3_part.R --task "$SLURM_ARRAY_TASK_ID" --cores "$SLURM_CPUS_PER_TASK" $EXTRA_ARGS
