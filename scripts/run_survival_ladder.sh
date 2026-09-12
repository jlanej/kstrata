#!/usr/bin/env bash
# k ladder over the spotlight regions (small query) against the seven genomes: input to the survival-length track.
set -euo pipefail
cd "$(dirname "$0")/.."
K=target/release/kstrata; P=out/prep; R=out/runs; mkdir -p $R
SUBJ="--subject hg002=$P/hg002 --subject chimp=$P/chimp --subject bonobo=$P/bonobo --subject gorilla=$P/gorilla --subject borang=$P/borang --subject sorang=$P/sorang --subject siamang=$P/siamang"
for k in 16 21 25 31 41 51 61 81 101 151 201 301 501 701 1001 1501 2001; do
  [ -s $R/regions_k${k}.json ] && continue
  $K run --query $P/regions --query-name regions $SUBJ -k $k --stride 1 --small-query --no-self --out $R/regions_k${k} 2>&1 | grep -E 'done in'
  echo "=== k$k done $(date)"
done
echo "=== LADDER DONE"
