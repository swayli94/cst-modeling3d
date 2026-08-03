#!/bin/bash
# Run every example script in its own directory against the repo package.
#
# Usage:  ./scripts/run_all_examples.sh [target_dir]
#
# Without arguments, runs the repo's own example/ tree IN PLACE (generated
# outputs are git-ignored). With a target_dir, expects an example tree copy
# there and runs it instead. Results are written to <target>/_results.txt.
#
# Environment: uses the repo's cst_modeling via PYTHONPATH and a headless
# matplotlib backend.
REPO=$(cd "$(dirname "$0")/.." && pwd)
ROOT=${1:-$REPO/example}
export PYTHONPATH=$REPO
export MPLBACKEND=Agg
LOG=$ROOT/_results.txt
: > $LOG
find $ROOT -name '*.py' | sort | while read f; do
  d=$(dirname $f); b=$(basename $f)
  out=$(cd $d && timeout 180 python $b 2>&1)
  code=$?
  echo "=== ${f#$ROOT/} | exit=$code" >> $LOG
  if [ $code -ne 0 ]; then
    echo "$out" | tail -5 | sed 's/^/    /' >> $LOG
  fi
done
echo "results written to $LOG"
