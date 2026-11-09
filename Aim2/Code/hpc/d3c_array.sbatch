#!/bin/bash
# =============================================================================
# hpc/d3c_array.sbatch
# Purpose : D3c on the cluster, one array task per cohort x variant
#           (10 cohorts x 2 variants = 20 tasks with the full settings).
# Usage   : from the project folder that holds Code/ (and Rlib/ with Rcpp):
#             mkdir -p logs && sbatch Code/hpc/d3c_array.sbatch
#           Then copy Results/test/objects/D3c_parts/ back and run
#             Rscript Code/run_all.R --steps D3c
# Container: R 4.4.2 image; Rcpp is not in the image and is loaded from
#           ./Rlib (installed once with the image's R).
# =============================================================================
#SBATCH -J aim2_d3c
#SBATCH -p red
#SBATCH --qos=xxlarge
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -a 1-20
#SBATCH -o logs/d3c_%A_%a.out

IMG=${IMG:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif}
export SINGULARITYENV_R_LIBS="$PWD/Rlib"
export APPTAINERENV_R_LIBS="$PWD/Rlib"
singularity exec --cleanenv -B /share,/tmp "$IMG" \
  Rscript Code/hpc/d3c_part.R --task "$SLURM_ARRAY_TASK_ID" --cores "$SLURM_CPUS_PER_TASK" $EXTRA_ARGS
