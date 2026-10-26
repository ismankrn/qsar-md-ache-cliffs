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
# A wedged link hangs scp forever, which stalled staging on a 1 KB file. These make
# a dead connection fail in seconds instead.
SSHOPT=(-o BatchMode=yes -o ConnectTimeout=15 -o ServerAliveInterval=15
        -o ServerAliveCountMax=3)
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MANIFEST=${MANIFEST:-$ROOT/results/14_md_runs_phase1.csv}
# BRIN is no longer in the chain: grompp runs on the instance, so nothing is
# pulled from the cluster. The variable is kept out deliberately rather than
# left dangling.
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

  ssh -n "${SSHOPT[@]}" "$H" 'mkdir -p /root/qsar/{tpr,systems,mdp,work,out,reduced,logs}'
  # NSTEPS_PROD rides along deliberately. Without it the worker finds no override, drops
  # -nsteps, and runs the tpr's full 25,000,000 steps: 50 ns instead of the 25 ns the rest
  # of the project used, which would look fine and not be comparable with anything.
  scp -q "${SSHOPT[@]}" "$ROOT/scripts/rung4_descriptors.py" "$ROOT/md/vast_worker.sh" \
         "$ROOT/md/NSTEPS_PROD" "$ROOT/docking/4ey7_rec_md.pdb" "$SLICE" "$H:/root/qsar/"
  ssh -n "${SSHOPT[@]}" "$H" "mv /root/qsar/$(basename $SLICE) /root/qsar/slice.txt"
  scp -q "${SSHOPT[@]}" "$ROOT"/md/mdp/{em,npt,prod}.mdp "$H:/root/qsar/mdp/"

  # Results already in hand go BACK to the instance. The worker skips any run whose JSON
  # exists, and that is the only thing stopping a freshly staged instance from recomputing
  # trajectories that have already been paid for. It matters whenever instances are swapped
  # mid-campaign: on 2026-10-23 three 250 GB instances were replaced by 32 GB ones, and the
  # slice assignment is a stride over the manifest, so the same runs land on the same slice
  # index again -- without this they would simply have been run a second time.
  back=0
  while read -r R; do
    J="$ROOT/results/16_rung4_per_run/$R.json"
    [ -f "$J" ] || continue
    scp -q "${SSHOPT[@]}" "$J" "$H:/root/qsar/out/" && back=$((back+1))
  done < "$SLICE"
  [ "$back" -gt 0 ] && echo "   JSON lama dikirim balik: $back run akan dilewati"

  # Systems rather than tprs. grompp happens on the instance, so nothing has to exist on
  # a third machine first; see vast_worker.sh. Only this instance's share is sent.
  miss=0
  while read -r R; do
    # ssh -n: without it ssh reads the loop's stdin and swallows the rest of the slice
    # file, so only the first run is ever staged.
    CID=${R%_r*}
    if ssh -n "${SSHOPT[@]}" "$H" "[ -f /root/qsar/tpr/$R.prod.tpr ] || [ -f /root/qsar/systems/$CID/system.gro ]" 2>/dev/null; then
      :
    else
      ssh -n "${SSHOPT[@]}" "$H" "mkdir -p /root/qsar/systems/$CID" \
        && scp -q "${SSHOPT[@]}" "$ROOT/md/systems/$CID/system.gro" "$ROOT/md/systems/$CID/system.top" \
              "$H:/root/qsar/systems/$CID/" \
        || { echo "   ! sistem gagal: $CID"; miss=$((miss+1)); continue; }
    fi
    # the per-run velocity seed lives in nvt.mdp, so it is sent per run and not per ligand
    scp -q "${SSHOPT[@]}" "$ROOT/md/runs/$R/nvt.mdp" "$H:/root/qsar/mdp/$R.nvt.mdp" 2>/dev/null \
      || { echo "   ! nvt.mdp gagal: $R"; miss=$((miss+1)); }
  done < "$SLICE"
  have=$(ssh -n "${SSHOPT[@]}" "$H" 'ls -d /root/qsar/systems/*/ 2>/dev/null | wc -l')
  echo "   sistem di instance: $have / $n   (gagal $miss)"
done
echo
echo "siap. jalankan di tiap instance:  ssh <host> 'cd /root/qsar && nohup bash vast_worker.sh > logs/worker.log 2>&1 &'"
