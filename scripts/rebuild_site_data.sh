#!/usr/bin/env bash
# After scripts/fetch_data.sh: prepare genomes, run k=31 at base resolution and k=1001 on the 1/8 sample, export the site data.
set -euo pipefail
cd "$(dirname "$0")/.."
K=target/release/kstrata; P=out/prep; R=out/runs; mkdir -p $P $R
until grep -q ALL_DATA_DONE data/fetch.log 2>/dev/null; do sleep 20; done
echo "=== data ready $(date)"
[ -s $P/chm13.idx ] || $K prep data/chm13v2.0.fa.gz $P/chm13
[ -s $P/hg002.idx ] || $K prep data/hg002v1.1.fasta.gz $P/hg002
for a in chimp bonobo gorilla borang sorang siamang; do [ -s $P/$a.idx ] || $K prep data/$a.*.fna.gz $P/$a; done
echo "=== prepped $(date)"
SUBJ="--subject hg002=$P/hg002 --subject chimp=$P/chimp --subject bonobo=$P/bonobo --subject gorilla=$P/gorilla --subject borang=$P/borang --subject sorang=$P/sorang --subject siamang=$P/siamang"
[ -s $R/chm13_k31_all.json ] || /usr/bin/time -l $K run --query $P/chm13 --query-name chm13 $SUBJ -k 31 --stride 1 --budget-gb 4 --out $R/chm13_k31_all 2>&1 | grep -E 'done in|peak memory'
echo "=== k31 done $(date)"
[ -s $R/chm13_k1001_all.json ] || /usr/bin/time -l $K run --query $P/chm13 --query-name chm13 $SUBJ -k 1001 --stride 8 --budget-gb 5 --out $R/chm13_k1001_all 2>&1 | grep -E 'done in|peak memory'
echo "=== k1001 done $(date)"
.venv/bin/python scripts/export_site.py 2>&1 | tail -3
echo "=== SITE DATA DONE $(date)"; du -sh docs/data
