#!/bin/bash
# =============================================================================
# hpc/m1_chain.sh RUNS - after a run is assembled: decide the next Model1 runs
# (hpc/m1_chain.R; wave 2, V1 folds + phenotype-term sensitivity, S2/S3) and launch them
# with hpc/de_launch.sh. A mkdir lock keeps two finishing runs from launching twice.
# Run from the task folder.
# =============================================================================
RUNS="$1"
source Code/hpc/de_common.sh
mkdir -p "$RUNS/_chain"
LOCK="$RUNS/_chain/lock"
for i in $(seq 1 120); do mkdir "$LOCK" 2>/dev/null && break; sleep 5; done
singularity exec --cleanenv -B /share,/tmp "$IMG" Rscript Code/hpc/m1_chain.R --runs "$RUNS"
TODO="$RUNS/_chain/todo.txt"
if [ -s "$TODO" ]; then
  while read -r RUN MODEL DATA WP WM W S FOLD; do
    nz() { [ "$1" = "-" ] && echo "" || echo "$1"; }
    WARM_POP=$(nz "$WP") WARM_MAPS=$(nz "$WM") PHENO_W=$(nz "$W") PHENO_S=$(nz "$S") FOLD=$(nz "$FOLD") \
      bash Code/hpc/de_launch.sh "$RUN" "$MODEL" "$DATA" >> "$RUNS/_chain/launched.log"
    echo "$(date '+%F %T') $RUN $MODEL data=$DATA warm=$WP w=$W s=$S fold=$FOLD" >> "$RUNS/_chain/launched.log"
  done < "$TODO"
  rm -f "$TODO"
fi
rmdir "$LOCK"
