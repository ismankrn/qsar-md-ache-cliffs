#!/bin/bash
#SBATCH --job-name=md_extract
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
# Rung-4 extraction for runs whose production finished on BRIN.
#
# This cannot run on the login node: the first extraction succeeded in 90 s and the second
# was SIGKILLed by the node's watchdog. Ask for what is needed and nothing more -- 1 node and
# 8 cores backfills, where a large request waits behind the queue.
#
# Four extractions at a time: each is single-threaded and reads one 2 GB trajectory, so the
# limit is I/O rather than CPU, and 8 cores leaves room for the trjconv calls.
set -u
module load bioinformatics/gromacs/2023-plumed-mpi 2>/dev/null || \
module load bioinformatics/gromacs/2023.3-mpi
GMX=$(command -v gmx 2>/dev/null || command -v gmx_mpi 2>/dev/null)
# The GROMACS module sets PYTHONHOME to the cluster's own Python, which breaks conda's
# interpreter outright ("Failed to import encodings module"). Drop it before activating.
unset PYTHONHOME PYTHONPATH
source ~/miniforge3/etc/profile.d/conda.sh
conda activate qsar-md
python -c "import numpy, pandas, MDAnalysis" || { echo "FATAL: env python tidak lengkap"; exit 1; }

ROOT=$HOME/qsar-md
RUNS=$ROOT/md/runs
OUT=$ROOT/results/16_rung4_per_run
RED=$ROOT/md/reduced
mkdir -p "$OUT" "$RED"
NPAR=${NPAR:-4}

one () {
  local d=$1 R
  R=$(basename "$d")
  [ -f "$OUT/$R.json" ] && return 0
  python "$ROOT/scripts/rung4_descriptors.py" --only "$R" --runs-dir "$RUNS" \
      --receptor "$ROOT/docking/4ey7_rec_md.pdb" --out "$OUT" > "$d/extract.log" 2>&1 || {
    echo "  $R: EKSTRAKSI GAGAL"; tail -3 "$d/extract.log" | sed 's/^/     /'; return 1; }
  local TPR="$d/prod.tpr"
  printf '"Protein" | "UNK"\nq\n' | "$GMX" make_ndx -f "$TPR" -o "$d/sel.ndx" > /dev/null 2>&1
  local GRP
  GRP=$(grep -o '^\[ [^]]* \]' "$d/sel.ndx" | tail -1 | tr -d '[] ')
  if [ -n "$GRP" ]; then
    "$GMX" trjconv -s "$TPR" -f "$d/prod.xtc" -o "$RED/$R.xtc" -dt 100 \
        -n "$d/sel.ndx" -pbc mol <<< "$GRP" > /dev/null 2>&1
    "$GMX" trjconv -s "$TPR" -f "$d/prod.xtc" -o "$RED/$R.gro" -dump 0 \
        -n "$d/sel.ndx" -pbc mol <<< "$GRP" > /dev/null 2>&1
  fi
  [ -f "$d/prod.cpt" ] && cp "$d/prod.cpt" "$RED/$R.cpt"
  echo "  $R: OK"
}

export -f one; export OUT RED RUNS ROOT GMX
n=0
# Gate on the step the run actually reached, read from prod.log, not on the trajectory's size.
# The old floor was MIN_XTC_MB=400, written when BRIN's tpr files were 50 ns and the in-flight
# runs were already past 25 ns, so "400 MB" then meant "well beyond the length we need". The
# 25 ns decision inverted that reasoning without anyone editing the number: a complete 25 ns
# prod.xtc here is about 890 MB (measured on job 555988, 146 MB per 4.1 ns), so 400 MB now
# means "accept at 11 ns". Extracting early is worse than skipping, because every later pass
# skips a run that already has a JSON -- a half-length result would become permanent, and
# silent except for a QC flag that is easy to read as a simulation problem.
MIN_STEPS=${MIN_STEPS:-$(cat "$ROOT/md/NSTEPS_PROD" 2>/dev/null || echo 12500000)}
MIN_STEPS=$(echo "$MIN_STEPS" | tr -dc '0-9')
echo "gerbang: step >= $MIN_STEPS, atau prod.gro ada"
for d in "$RUNS"/*/; do
  [ -f "$d/prod.xtc" ] || continue
  R=$(basename "$d"); [ -f "$OUT/$R.json" ] && continue
  # Two conditions, and both are needed.
  #
  # Far enough: the last "Step Time" block in the log must be at or past the floor. A run cut
  # short by -maxh below the floor is skipped, which is what we want -- it is resumable with
  # -cpi and worth resuming rather than analysing short.
  st=$(grep -A1 -E '^ *Step *Time' "$d/prod.log" 2>/dev/null | grep -oE '^ *[0-9]+' | tail -1 | tr -d ' ')
  [ -z "$st" ] && continue
  [ "$st" -lt "$MIN_STEPS" ] && { echo "  $R: $st step, belum $MIN_STEPS -- dilewati"; continue; }
  #
  # And finished writing: mdrun has either written its final structure or logged that it is
  # done. The step count alone is not enough. A run that has just reached the target is still
  # flushing its trajectory, and reading it then yields a short file that looks like a real
  # result -- the same silent, permanent failure the step floor was introduced to prevent,
  # arriving by a different route. Measured on job 555988: at 12,450,000 of 12,500,000 steps a
  # run had neither marker; at 12,500,000 it had both, with prod.xtc at 846 MB.
  #
  # prod.gro and the log line are not redundant. mdrun writes confout only when it runs the
  # full nsteps, so a -maxh stop leaves no prod.gro at all -- all sixteen of the 50 ns runs
  # cut at the walltime boundary lack one, and all sixteen logged "Finished mdrun". Testing
  # only for prod.gro would silently drop every one of them.
  if [ ! -f "$d/prod.gro" ] && ! grep -q "Finished mdrun" "$d/prod.log" 2>/dev/null; then
    echo "  $R: $st step tapi mdrun belum tuntas menulis -- dilewati"
    continue
  fi
  one "$d" &
  n=$((n+1))
  [ $((n % NPAR)) -eq 0 ] && wait
done
wait
echo "--- JSON: $(ls "$OUT"/*.json 2>/dev/null | wc -l), arsip: $(ls "$RED"/*.xtc 2>/dev/null | wc -l) ---"
