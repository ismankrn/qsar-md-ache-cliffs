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
# Gate on the trajectory, not on prod.gro. A run stopped by -maxh at the walltime boundary
# exits cleanly but there is no guarantee of a confout, and BRIN's in-flight runs are already
# past 25 ns -- long enough to analyse and too expensive to discard over a missing file. The
# extraction's own min-frames check is the real floor; this only avoids opening empty runs.
MIN_XTC_MB=${MIN_XTC_MB:-400}
for d in "$RUNS"/*/; do
  [ -f "$d/prod.xtc" ] || continue
  if [ ! -f "$d/prod.gro" ]; then
    sz=$(du -m "$d/prod.xtc" 2>/dev/null | cut -f1)
    [ "${sz:-0}" -lt "$MIN_XTC_MB" ] && continue
  fi
  R=$(basename "$d"); [ -f "$OUT/$R.json" ] && continue
  one "$d" &
  n=$((n+1))
  [ $((n % NPAR)) -eq 0 ] && wait
done
wait
echo "--- JSON: $(ls "$OUT"/*.json 2>/dev/null | wc -l), arsip: $(ls "$RED"/*.xtc 2>/dev/null | wc -l) ---"
