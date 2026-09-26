#!/bin/bash
# Drive one batch of MD runs through em -> nvt -> npt -> prod.
#
# Resumable at every level, because a 50 ns run at ~35 ns/day takes 1.4 days and a job
# that dies at hour 30 must not start over:
#   - a stage whose .gro already exists is skipped
#   - production restarts from prod.cpt via -cpi, and -maxh makes mdrun checkpoint and
#     exit cleanly before the walltime kills it mid-write
#
# Runs are driven with `gmx_mpi mdrun -multidir`, one MPI rank per simulation and
# -ntomp threads each, so a single SLURM job carries many simulations. With MaxJobs 3
# on this cluster, few large jobs are the only way to get real parallelism.
#
# Usage: BATCH=<i> NBATCH=<n> NTOMP=<threads> bash md/run_md.sh
set -u
ROOT="$HOME/qsar-md"
RUNS="$ROOT/md/runs"
NTOMP=${NTOMP:-16}
BATCH=${BATCH:-0}
NBATCH=${NBATCH:-1}
MAXH=${MAXH:-68}          # leave headroom inside a 72 h walltime

# The MPI GROMACS module installs only `gmx_mpi`; the serial one only `gmx`. grompp runs
# single-rank either way, so resolve whichever exists instead of hard-coding one. Hard-coding
# `gmx` is what made all 379 grompp calls fail under the MPI module, leaving no .tpr for
# mdrun to read.
GROMPP=$(command -v gmx 2>/dev/null || command -v gmx_mpi 2>/dev/null)
if [ -z "$GROMPP" ]; then
  echo "FATAL: tidak ada gmx maupun gmx_mpi di PATH -- modul GROMACS belum dimuat?"
  exit 1
fi
echo "grompp memakai: $GROMPP"

# MANIFEST selects which runs this job covers. Phase 1 is replicate 1 only (217 runs);
# the 162 replicate trajectories are deferred until phase 1 shows a signal, because a null
# result needs no reproducibility check. Their systems are already equilibrated and waiting.
MANIFEST=${MANIFEST:-$ROOT/results/14_md_runs_phase1.csv}
if [ ! -f "$MANIFEST" ]; then echo "FATAL: manifest tidak ada: $MANIFEST"; exit 1; fi
echo "manifest: $(basename $MANIFEST) ($(( $(wc -l < "$MANIFEST") - 1 )) run)"
mapfile -t ALL < <(tail -n +2 "$MANIFEST" | cut -d, -f1)
DIRS=()
for i in $(seq $BATCH $NBATCH $((${#ALL[@]}-1))); do
  DIRS+=("$RUNS/${ALL[$i]}")
done
echo "batch $BATCH/$NBATCH -> ${#DIRS[@]} run, $NTOMP thread per simulasi"

# ---- skip runs that already finished production
TODO=()
for d in "${DIRS[@]}"; do
  [ -f "$d/prod.gro" ] && continue
  TODO+=("$d")
done
echo "belum selesai: ${#TODO[@]}"
[ ${#TODO[@]} -eq 0 ] && { echo "tidak ada yang dikerjakan"; exit 0; }

NRANK=${#TODO[@]}

# ---- grompp is serial and cheap; do it per directory, per stage, only if needed
prep () {   # stage  coord_from  extra
  local stage=$1 coord=$2; shift 2
  local nfail=0
  for d in "${TODO[@]}"; do
    [ -f "$d/$stage.tpr" ] && continue
    [ -f "$d/$stage.gro" ] && continue
    ( cd "$d" && "$GROMPP" grompp -f $stage.mdp -c "$coord" -p system.top \
        -o $stage.tpr -maxwarn 10 "$@" > grompp_$stage.log 2>&1 ) \
      || { echo "  GROMPP GAGAL $stage di $(basename $d)"; tail -5 "$d/grompp_$stage.log"; nfail=$((nfail+1)); }
  done
  # Fail the stage rather than handing mdrun a directory with no .tpr. That is what turned
  # one missing binary into 379 identical "File input/output error: em.tpr" messages.
  if [ $nfail -gt 0 ]; then
    echo "  $stage: $nfail grompp gagal -- berhenti sebelum mdrun"
    return 1
  fi
}

# Never launch more ranks than SLURM allocated: -multidir takes one rank per simulation,
# so a batch holding more runs than there are task slots would silently oversubscribe the
# nodes and every simulation in it would crawl. Process in chunks of at most SLURM_NTASKS.
MAXSIM=${SLURM_NTASKS:-${#TODO[@]}}

run () {    # stage  [extra mdrun args]
  local stage=$1; shift
  local pend=()
  for d in "${TODO[@]}"; do
    [ -f "$d/$stage.gro" ] || pend+=("$d")
  done
  [ ${#pend[@]} -eq 0 ] && { echo "  $stage: sudah lengkap"; return 0; }
  echo "  $stage: ${#pend[@]} simulasi, maksimum $MAXSIM serentak"
  local i=0 chunk=0
  while [ $i -lt ${#pend[@]} ]; do
    local slice=("${pend[@]:$i:$MAXSIM}")
    chunk=$((chunk+1))
    echo "    potongan $chunk: ${#slice[@]} simulasi"
    # -maxh must reflect the time this job has LEFT, not a fixed number. Chunks run one
    # after another, so a later chunk starting with 20 h left would otherwise be told it
    # has 70 h, run past the walltime and be killed mid-stride -- losing up to a
    # checkpoint interval instead of stopping cleanly. SLURM reports the remainder.
    local left
    left=$(squeue -h -j "${SLURM_JOB_ID:-0}" -o %L 2>/dev/null)
    local h=$MAXH
    if [ -n "$left" ]; then
      h=$(awk -v t="$left" 'BEGIN{
            n=split(t,a,"-"); d=0; rest=t
            if(n==2){d=a[1]; rest=a[2]}
            split(rest,b,":")
            if(length(b)==3) s=b[1]*3600+b[2]*60+b[3]
            else if(length(b)==2) s=b[1]*60+b[2]
            else s=b[1]*60
            printf "%.2f", d*24 + s/3600 - 0.25}')
    fi
    awk -v h="$h" 'BEGIN{if(h<=0.1) exit 1}' || { echo "    waktu job habis, berhenti"; return 0; }
    echo "      -maxh $h"
    mpirun -np ${#slice[@]} gmx_mpi mdrun -multidir "${slice[@]}" \
        -deffnm $stage -ntomp $NTOMP -pin on -maxh $h "$@" \
        > "$ROOT/logs/mdrun_${stage}_b${BATCH}_c${chunk}.log" 2>&1
    local rc=$?
    [ $rc -ne 0 ] && echo "    potongan $chunk: mdrun keluar kode $rc (cek log)"
    i=$((i+MAXSIM))
  done
  return 0
}

echo "== minimisasi =="
prep em  system.gro || exit 1
run  em

echo "== NVT (restraint posisi) =="
# -r supplies the reference coordinates for POSRES; without it grompp fails, and with
# the wrong file the protein would be restrained to the wrong positions
prep nvt em.gro -r em.gro || exit 1
run  nvt

echo "== NPT =="
prep npt nvt.gro -r nvt.gro -t nvt.cpt || exit 1
run  npt

if [ "${STOP_AFTER:-prod}" = "npt" ]; then
  echo "== berhenti setelah NPT (STOP_AFTER=npt) =="
else
  echo "== produksi 50 ns =="
  prep prod npt.gro -t npt.cpt || exit 1
  run  prod -cpi prod.cpt
fi

echo
echo "== status akhir batch $BATCH =="
if [ "${STOP_AFTER:-prod}" = "npt" ]; then
  ok=0; bad=0
  for d in "${TODO[@]}"; do
    if [ -f "$d/npt.gro" ]; then ok=$((ok+1)); else bad=$((bad+1)); echo "  $(basename $d): ekuilibrasi TIDAK selesai"; fi
  done
  echo "  ekuilibrasi selesai: $ok, gagal: $bad"
  exit 0
fi
for d in "${TODO[@]}"; do
  n=$(basename "$d")
  if [ -f "$d/prod.gro" ]; then echo "  $n: SELESAI"
  elif [ -f "$d/prod.cpt" ]; then
    t=$(gmx check -f "$d/prod.cpt" 2>&1 | sed -n 's/.*Last frame.*time *\([0-9.]*\).*/\1/p' | tail -1)
    echo "  $n: berjalan, checkpoint di ${t:-?} ps"
  else echo "  $n: belum sampai produksi"
  fi
done
