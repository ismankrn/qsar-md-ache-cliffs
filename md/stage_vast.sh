#!/bin/bash
# Split the run manifest across N rented instances and stage each one's share.
#
# Slices are disjoint by construction (stride), so two instances can never pick up the same
# run. That matters because the instances share no filesystem: there is no lock to fall back
# on, and a duplicated trajectory is wasted GPU-hours while a missed one is a hole in the
# dataset that only shows up at analysis time.
#
# Usage: bash md/stage_vast.sh <host1> [host2 ...]
#   hosts are ~/.ssh/config aliases; the manifest is split evenly across them.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MANIFEST=${MANIFEST:-$ROOT/results/14_md_runs_phase1.csv}
BRIN=${BRIN:-mahameru}
HOSTS=("$@")
[ ${#HOSTS[@]} -eq 0 ] && { echo "pakai: bash md/stage_vast.sh vast1 [vast2 ...]"; exit 1; }
# TOTAL fixes how many slices the manifest is cut into, independent of how many hosts exist
# right now. Staging one instance as slice 0 of 8 means the seven that arrive later take
# slices 1-7 and nothing is ever run twice -- the instances share no filesystem and there is
# no lock to fall back on.
N=${TOTAL:-${#HOSTS[@]}}
[ ${#HOSTS[@]} -gt $N ] && { echo "FATAL: host (${#HOSTS[@]}) melebihi TOTAL ($N)"; exit 1; }

# No mapfile: macOS ships bash 3.2, where it does not exist. This script runs on the
# laptop, unlike run_md.sh which runs on Linux and can use it.
RUNS=()
while IFS= read -r line; do RUNS+=("$line"); done < <(tail -n +2 "$MANIFEST" | cut -d, -f1)
echo "manifest : $(basename $MANIFEST), ${#RUNS[@]} run"
echo "instance : $N  -> ~$(( (${#RUNS[@]} + N - 1) / N )) run masing-masing"
echo

STAGE=$(mktemp -d); trap 'rm -rf "$STAGE"' EXIT
for i in "${!HOSTS[@]}"; do
  H="${HOSTS[$i]}"
  SLICE="$STAGE/slice_$i.txt"; : > "$SLICE"
  for ((j=i; j<${#RUNS[@]}; j+=N)); do echo "${RUNS[$j]}" >> "$SLICE"; done
  n=$(wc -l < "$SLICE")
  echo "== $H : $n run =="

  ssh -n -o BatchMode=yes "$H" 'mkdir -p /root/qsar/{tpr,work,out,reduced,logs}' 
  scp -q "$ROOT/scripts/rung4_descriptors.py" "$ROOT/md/vast_worker.sh" \
         "$ROOT/docking/4ey7_rec_md.pdb" "$SLICE" "$H:/root/qsar/" 
  ssh -n -o BatchMode=yes "$H" "mv /root/qsar/$(basename $SLICE) /root/qsar/slice.txt"

  # tpr files come from BRIN; pull only this instance's share rather than all 217
  miss=0
  while read -r R; do
    # ssh -n: without it ssh reads the loop's stdin and swallows the rest of the slice
    # file, so only the first run is ever staged.
    if ! ssh -n -o BatchMode=yes "$H" "[ -f /root/qsar/tpr/$R.prod.tpr ]" 2>/dev/null; then
      scp -q "$BRIN:~/qsar-md/md/runs/$R/prod.tpr" "$STAGE/$R.prod.tpr" 2>/dev/null \
        && scp -q "$STAGE/$R.prod.tpr" "$H:/root/qsar/tpr/$R.prod.tpr" \
        && rm -f "$STAGE/$R.prod.tpr" \
        || { echo "   ! tpr gagal: $R"; miss=$((miss+1)); }
    fi
  done < "$SLICE"
  have=$(ssh -n -o BatchMode=yes "$H" 'ls /root/qsar/tpr/*.prod.tpr 2>/dev/null | wc -l')
  echo "   tpr di instance: $have / $n   (gagal $miss)"
done
echo
echo "siap. jalankan di tiap instance:  ssh <host> 'cd /root/qsar && nohup bash vast_worker.sh > logs/worker.log 2>&1 &'"
