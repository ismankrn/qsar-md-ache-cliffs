#!/bin/bash
# Pull finished results off the rented instances.
#
# Everything a worker produces lives on its instance disk and nowhere else. Vast instances
# can be reclaimed, run out of rental time, or simply fail, and at ~6.3 h per trajectory an
# instance that dies late can take a week of GPU with it. Run this often; it is cheap
# because only descriptors and reduced archives move, never the full trajectories.
#
# Usage: bash md/pull_vast.sh vast1 [vast2 ...]
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST_JSON="$ROOT/results/16_rung4_per_run"
DEST_RED="$ROOT/md/reduced"
mkdir -p "$DEST_JSON" "$DEST_RED"
[ $# -eq 0 ] && { echo "pakai: bash md/pull_vast.sh vast1 [vast2 ...]"; exit 1; }

total_new=0
for H in "$@"; do
  n_remote=$(ssh -n -o BatchMode=yes -o ConnectTimeout=20 "$H" \
    'ls /root/qsar/out/*.json 2>/dev/null | wc -l' 2>/dev/null)
  if [ -z "$n_remote" ]; then echo "  $H: tidak terjangkau"; continue; fi
  before=$(ls "$DEST_JSON"/*.json 2>/dev/null | wc -l | tr -d ' ')
  # --update, NOT --ignore-existing. A run that has to be redone -- because its descriptors
  # were extracted by a defective version of the script, which has already happened once --
  # produces a JSON of the same name that must replace the stale one. --ignore-existing would
  # silently keep the bad file forever. --update still refuses to go backwards in time, so a
  # corrected local file is never clobbered by an older remote one.
  moved_j=$(rsync -a --update -i -e ssh "$H:/root/qsar/out/" "$DEST_JSON/" 2>/dev/null \
            | grep -c '^[<>]')
  rsync -a --update -e ssh "$H:/root/qsar/reduced/" "$DEST_RED/" 2>/dev/null
  after=$(ls "$DEST_JSON"/*.json 2>/dev/null | wc -l | tr -d ' ')
  fresh=$((after-before)); upd=$((moved_j-fresh)); [ $upd -lt 0 ] && upd=0
  slice=$(ssh -n -o BatchMode=yes "$H" 'wc -l < /root/qsar/slice.txt' 2>/dev/null)
  printf "  %-8s instance %s/%s selesai, +%s baru, %s diperbarui\n" \
    "$H:" "$n_remote" "${slice:-?}" "$fresh" "$upd"
  total_new=$((total_new + fresh))
done
echo "  ---"
echo "  JSON di laptop : $(ls "$DEST_JSON"/*.json 2>/dev/null | wc -l | tr -d ' ') / 217"
echo "  arsip          : $(ls "$DEST_RED"/*.xtc 2>/dev/null | wc -l | tr -d ' ')  ($(du -sh "$DEST_RED" 2>/dev/null | cut -f1))"
echo "  baru kali ini  : $total_new"
