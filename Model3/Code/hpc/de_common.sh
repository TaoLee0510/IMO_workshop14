# hpc/de_common.sh - sourced by the Model3 DE-ITS sbatch scripts. Task folder = $PWD.
IMG=${IMG:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif}
RLIB=${RLIB:-$PWD/../Rlib}
export SINGULARITYENV_R_LIBS="$RLIB"
export APPTAINERENV_R_LIBS="$RLIB"
rtask() { singularity exec --cleanenv -B /share,/tmp "$IMG" Rscript Code/hpc/de_task.R "$@"; }
# nodes where singularity could not resolve the user id ("unknown userid", 2026-11-05): excluded from every job
export SBATCH_EXCLUDE=hpctpa3pc0003,hpctpa3pc0005,hpctpa3pc0007,hpctpa3pc0011,hpctpa3pc0029,hpctpa3pc0051,hpctpa3pc0065
# --qos must be on the command line: the #SBATCH --qos line is not honoured on this cluster
sub() { sbatch --parsable --qos=xxlarge -p red "$@"; }
