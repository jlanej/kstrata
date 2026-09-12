#!/usr/bin/env bash
# Flagship k=31 run at base resolution, then the k ladder on a 1/8 systematic sample, then summaries.
set -euo pipefail
cd "$(dirname "$0")/.."
BIN=target/release/kstrata; PREP=out/prep; RUNS=out/runs; PY=.venv/bin/python
SUBJ="--subject hg002=$PREP/hg002 --subject chimp=$PREP/chimp --subject bonobo=$PREP/bonobo --subject gorilla=$PREP/gorilla --subject borang=$PREP/borang --subject sorang=$PREP/sorang --subject siamang=$PREP/siamang"
# wait for the six ape preps
until [ -s $PREP/chimp.idx ] && [ -s $PREP/bonobo.idx ] && [ -s $PREP/gorilla.idx ] && [ -s $PREP/borang.idx ] && [ -s $PREP/sorang.idx ] && [ -s $PREP/siamang.idx ] && ! pgrep -f 'kstrata prep' >/dev/null; do sleep 20; done
echo "=== apes ready $(date)"
if [ ! -s $RUNS/chm13_k31_all.json ]; then
  /usr/bin/time -l $BIN run --query $PREP/chm13 --query-name chm13 $SUBJ -k 31 --stride 1 --budget-gb 4 --out $RUNS/chm13_k31_all 2>&1 | grep -E 'k=31|subject siamang|done in|peak memory'
fi
echo "=== k31 done $(date)"
$PY scripts/strata.py $RUNS/chm13_k31_all --genome chm13 --window 10000 2>/dev/null | tail -3
$PY scripts/strata.py $RUNS/chm13_k31_all --genome chm13 --window 100000 2>/dev/null | tail -1
echo "=== k31 summaries done $(date)"
for k in 16 21 41 61 101 201 501 1001 2001; do
  if [ ! -s $RUNS/chm13_k${k}_all.json ]; then
    /usr/bin/time -l $BIN run --query $PREP/chm13 --query-name chm13 $SUBJ --subject random=$PREP/random -k $k --stride 8 --budget-gb 5 --out $RUNS/chm13_k${k}_all 2>&1 | grep -E "^\[.*k=|done in|peak memory"
  fi
  $PY scripts/strata.py $RUNS/chm13_k${k}_all --genome chm13 --window 10000 2>/dev/null | tail -1
  echo "=== k$k done $(date)"
done
echo "=== ALL DONE $(date)"
