#!/usr/bin/env bash
# Fetch every input: T2T-CHM13v2.0 + annotation, HG002 v1.1 + cenSat, and the six T2T ape assemblies (~8.4 GB).
# Idempotent: resumes partial files and retries until each file matches the server's size.
set -u
DATA="${1:-data}"; mkdir -p "$DATA"
HPP=https://s3-us-west-2.amazonaws.com/human-pangenomics
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA
get() {  # url dest
  local url=$1 dest=$2 want have
  want=$(curl -sIL --max-time 30 "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:" {n=$2} END {print n}')
  for attempt in $(seq 1 40); do
    have=$( [ -f "$dest" ] && stat -f %z "$dest" || echo 0 )
    if [ -n "$want" ] && [ "$have" -eq "$want" ]; then echo "  ok      $(basename "$dest")"; return 0; fi
    echo "  fetch   $(basename "$dest") (attempt $attempt, have $have of ${want:-?})"
    curl -fL --no-progress-meter --retry 5 --retry-delay 5 -C - -o "$dest" "$url" && [ -z "$want" ] && return 0
    sleep 5
  done
  echo "  FAILED  $url" >&2; return 1
}
echo "[1/4] T2T-CHM13v2.0 + annotation"
get $HPP/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz $DATA/chm13v2.0.fa.gz
for b in censat_v2.1 SD; do get $HPP/T2T/CHM13/assemblies/annotation/chm13v2.0_${b}.bed $DATA/chm13v2.0_${b}.bed; done
echo "[2/4] HG002 v1.1 + cenSat"
get $HPP/T2T/HG002/assemblies/hg002v1.1.fasta.gz $DATA/hg002v1.1.fasta.gz
get $HPP/T2T/HG002/assemblies/annotation/centromere/hg002v1.1_v2.0/hg002v1.1.cenSatv2.0.bed $DATA/hg002v1.1.cenSatv2.0.bed
echo "[3/4] six ape assemblies (Yoo et al. 2025, primary), three at a time"
ape() { get "$NCBI/$2/$(basename "$2")_genomic.fna.gz" "$DATA/$1.$(basename "$2").fna.gz"; }
ape chimp   028/858/775/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri &
ape bonobo  029/289/425/GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri &
ape gorilla 029/281/585/GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri &
wait
ape sorang  028/885/655/GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri &
ape borang  028/885/625/GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri &
ape siamang 028/878/055/GCA_028878055.2_NHGRI_mSymSyn1-v2.0_pri &
wait
echo "[4/4] done"; ls -lh "$DATA" | awk 'NR>1 {printf "  %-52s %8s\n", $9, $5}'; echo ALL_DATA_DONE
