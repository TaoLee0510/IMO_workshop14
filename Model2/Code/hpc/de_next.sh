#!/bin/bash
# =============================================================================
# hpc/de_next.sh RUN - read RUN/status and submit what comes next:
#   round0 0              -> E array of round 0 (independent fits) + round-0 initialisation (afterany)
#   missing0 0 i,j        -> rerun the missing round-0 tasks once, then the initialisation again
#   estep K | continue K  -> E array (round K) + M step (afterany)
#   converged K | cap K   -> final E array (round K, importance sampling) + assemble
#   missing K i,j         -> rerun the missing E tasks once, then the M step again
#   missing_final i,j     -> rerun the missing final tasks once, then assemble again
# Every submission is logged to RUN/submissions.log.
# =============================================================================
RUN="$1"
source Code/hpc/de_common.sh
read -r ST K LIST < "$RUN/status"
N=$(cat "$RUN/n_patients")
NAME=$(basename "$RUN")
mkdir -p "$RUN/logs"
logsub() { echo "$(date '+%F %T') $*" >> "$RUN/submissions.log"; }
round() {   # $1 = WHAT, $2 = K, $3 = array spec
  local e
  e=$(sub --array="$3" -J "e_${NAME}" -o "$RUN/logs/${1}_${2}_%A_%a.out" \
      --export=ALL,RUN="$RUN",ITER="$2",WHAT="$1" Code/hpc/de_estep.sbatch)
  logsub "$1 round $2 array $3 job $e"
  echo "$e"
}
case "$ST" in
  round0|missing0)
    if [ "$ST" = "missing0" ]; then
      if [ -e "$RUN/retry_0" ]; then logsub "STOP: round 0 still missing $LIST"; echo "error 0 missing $LIST" > "$RUN/status"; exit 1; fi
      touch "$RUN/retry_0"; ARR="$LIST"
    else ARR="1-$N"; fi
    e=$(round estep 0 "$ARR")
    r=$(sub --dependency=afterany:"$e" -J "r0_${NAME}" -o "$RUN/logs/round0_%j.out" \
        --export=ALL,RUN="$RUN",WHAT=round0 Code/hpc/de_small.sbatch)
    logsub "round0 init job $r" ;;
  estep|continue)
    e=$(round estep "$K" "1-$N")
    m=$(sub --dependency=afterany:"$e" -J "m_${NAME}" -o "$RUN/logs/mstep_${K}_%j.out" \
        --export=ALL,RUN="$RUN",ITER="$K" Code/hpc/de_mstep.sbatch)
    logsub "mstep round $K job $m" ;;
  converged|cap)
    e=$(round efinal "$K" "1-$N")
    a=$(sub --dependency=afterany:"$e" -J "a_${NAME}" -o "$RUN/logs/assemble_%j.out" \
        --export=ALL,RUN="$RUN",WHAT=assemble Code/hpc/de_small.sbatch)
    logsub "assemble job $a ($ST)" ;;
  missing)
    if [ -e "$RUN/retry_$K" ]; then logsub "STOP: round $K still missing $LIST"; echo "error $K missing $LIST" > "$RUN/status"; exit 1; fi
    touch "$RUN/retry_$K"
    e=$(round estep "$K" "$LIST")
    m=$(sub --dependency=afterany:"$e" -J "m_${NAME}" -o "$RUN/logs/mstep_${K}_%j.out" \
        --export=ALL,RUN="$RUN",ITER="$K" Code/hpc/de_mstep.sbatch)
    logsub "retry mstep round $K job $m" ;;
  missing_final)
    if [ -e "$RUN/retry_final" ]; then logsub "STOP: final still missing $K"; echo "error final missing $K" > "$RUN/status"; exit 1; fi
    touch "$RUN/retry_final"
    KF=$(ls "$RUN"/pop_*.rds | sed 's/.*pop_0*//; s/\.rds//' | sort -n | tail -1)
    e=$(round efinal "$KF" "$K")
    a=$(sub --dependency=afterany:"$e" -J "a_${NAME}" -o "$RUN/logs/assemble_%j.out" \
        --export=ALL,RUN="$RUN",WHAT=assemble Code/hpc/de_small.sbatch)
    logsub "retry assemble job $a" ;;
  *) logsub "no action for status: $ST $K $LIST" ;;
esac
