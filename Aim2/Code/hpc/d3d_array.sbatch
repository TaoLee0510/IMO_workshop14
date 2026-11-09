#!/bin/bash
# =============================================================================
# hpc/d3d_array.sbatch
# Purpose : D3d on the cluster, one array task per cohort x variant
#           (50 cohorts x 2 variants = 100 tasks with the full settings).
# Usage   : from the project folder that holds Code/, Rlib/ and
#           Results/real/objects/D3d_truth.rds (no patient data needed):
#             mkdir -p logs && sbatch Code/hpc/d3d_array.sbatch
#           Then copy Results/real/objects/D3d_parts/ back and run
#             Rscript Code/run_all.R --data real --steps D3d
# Container: as d3c_array.sbatch (R 4.4.2 image, Rcpp from ./Rlib).
# =============================================================================
#SBATCH -J aim2_d3d
#SBATCH -p red
#SBATCH --qos=xxlarge
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -a 1-100
#SBATCH -o logs/d3d_%A_%a.out

IMG=${IMG:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif}
export SINGULARITYENV_R_LIBS="$PWD/Rlib"
export APPTAINERENV_R_LIBS="$PWD/Rlib"
singularity exec --cleanenv -B /share,/tmp "$IMG" \
  Rscript Code/hpc/d3d_part.R --task "$SLURM_ARRAY_TASK_ID" --cores "$SLURM_CPUS_PER_TASK" $EXTRA_ARGS
