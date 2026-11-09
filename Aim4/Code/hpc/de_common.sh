# hpc/de_common.sh - sourced by the DE-ITS sbatch scripts (v1.3). Task folder = $PWD.
IMG=${IMG:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif}
RLIB=${RLIB:-$PWD/../Rlib}
export SINGULARITYENV_R_LIBS="$RLIB"
export APPTAINERENV_R_LIBS="$RLIB"
rtask() { singularity exec --cleanenv -B /share,/tmp "$IMG" Rscript Code/hpc/de_task.R "$@"; }
# --qos must be on the command line: the #SBATCH --qos line is not honoured on this cluster
sub() { sbatch --parsable --qos=xxlarge -p red "$@"; }
