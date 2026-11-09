#!/bin/bash
# =============================================================================
# hpc/k3_array.sbatch
# Purpose : K3 on the cluster, one array task per synthetic cohort
#           (real35: 4 mechanisms x 2 clone settings x 20 replicates = 160;
#            n350: 5 replicates x 2 VAF variants = 10; total 170).
# Usage   : from the Aim4 task folder that holds Code/ and Data/derived/
#             mkdir -p logs && sbatch --qos=xxlarge --array=1-170 Code/hpc/k3_array.sbatch
#           (--qos must be given on the command line: the #SBATCH --qos line below is not
#            honoured on this cluster and the job lands in the 'small' QOS)
#           then copy Results/real/objects/K3_parts/ back and run
#             Rscript Code/run_all.R --data real --steps K3
# Container: R 4.4.2 image; Rcpp is loaded from ../Rlib (shared with Aim 2).
# =============================================================================
#SBATCH -J aim4_k3
#SBATCH -p red
#SBATCH --qos=xxlarge
#SBATCH -c 16
#SBATCH --mem=24G
#SBATCH -t 06:00:00
#SBATCH -o logs/k3_%A_%a.out

IMG=${IMG:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif}
RLIB=${RLIB:-$PWD/../Rlib}
export SINGULARITYENV_R_LIBS="$RLIB"
export APPTAINERENV_R_LIBS="$RLIB"
singularity exec --cleanenv -B /share,/tmp "$IMG" \
  Rscript Code/hpc/k3_part.R --task "$SLURM_ARRAY_TASK_ID" --cores "$SLURM_CPUS_PER_TASK" $EXTRA_ARGS
