#!/bin/bash
# One Vast instance working through its slice of the run manifest.
#
# Per run: mdrun 50 ns -> extract rung-4 descriptors -> keep the JSON and a reduced
# trajectory -> DELETE the full trajectory. That last step is what makes a 32 GB disk
# enough: the full xtc is ~2 GB per run and 34 runs would need 68 GB, but the descriptors
# are tens of kilobytes and the reduced archive ~50 MB. Without it every instance would
# need >100 GB and cost more.
#
# Resumable at two levels: a run whose JSON exists is skipped, and an interrupted mdrun
# restarts from prod.cpt via -cpi. Vast instances can be reclaimed, and at 6.3 h per
# trajectory an interruption must not cost the whole run.
#
# Usage: SLICE=<file with run names> bash vast_worker.sh
# No `set -u`: conda's GMXRC reads unset variables and aborts under it. This exact
# failure already cost a run on the campus GPU node; it is written down here so the
# next script does not rediscover it a third time.
ROOT=/root/qsar
SLICE=${SLICE:-$ROOT/slice.txt}
NTOMP=${NTOMP:-4}          # measured optimum: 4 threads beat 8 and 11 on this GPU
# mdrun refuses to start when OMP_NUM_THREADS is set to anything other than -ntomp:
# "Environment variable OMP_NUM_THREADS and the number of threads requested on the
# command line have different values." Some Vast images set it in /etc/environment to
# nproc, some leave it unset, so the same command that ran for days on one instance
# failed all 12 runs in 40 seconds on the next. Export it rather than hope.
export OMP_NUM_THREADS=$NTOMP
KEEP_REDUCED=${KEEP_REDUCED:-1}
# NSTEPS overrides the tpr's 25,000,000 (50 ns). Used only to prove the pipeline end to end
# on a short run: finding a bug in the extraction step after 6.3 h of mdrun is an expensive
# way to learn about it.
NSTEPS=${NSTEPS:-}
EXTRA_EXTRACT=${EXTRA_EXTRACT:-}

source /root/miniforge3/etc/profile.d/conda.sh
conda activate gmx
export GMX_MAXBACKUP=-1

mkdir -p "$ROOT"/{tpr,work,out,reduced,logs}
[ -f "$SLICE" ] || { echo "FATAL: slice tidak ada: $SLICE"; exit 1; }
echo "instance : $(hostname)   GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "slice    : $(wc -l < "$SLICE") run, $NTOMP thread"
echo "mulai    : $(date)"
echo

done_n=0; skip=0; fail=0; consec=0
while read -r RUN; do
  [ -z "$RUN" ] && continue
  if [ -f "$ROOT/out/$RUN.json" ]; then skip=$((skip+1)); continue; fi
  TPR="$ROOT/tpr/$RUN.prod.tpr"
  if [ ! -f "$TPR" ]; then echo "  $RUN: tpr tidak ada, dilewati"; fail=$((fail+1)); continue; fi

  W="$ROOT/work/$RUN"; mkdir -p "$W"
  # Link the tpr into the work directory. mdrun takes -s explicitly, but the extractor
  # resolves its topology as <run_dir>/prod.tpr and falls back to npt.gro, which does not
  # exist here -- that mismatch is what failed the first end-to-end attempt. Linking keeps
  # rung4_descriptors.py byte-identical between BRIN and the rented instances.
  ln -sf "$TPR" "$W/prod.tpr"
  t0=$SECONDS
  NS_ARG=""; [ -n "$NSTEPS" ] && NS_ARG="-nsteps $NSTEPS"
  ( cd "$W" && gmx mdrun -s "$TPR" -deffnm prod -ntmpi 1 -ntomp "$NTOMP" \
      -nb gpu -pme gpu -bonded gpu -update gpu -cpi prod.cpt -noconfout $NS_ARG \
      > mdrun.log 2>&1 )
  if [ ! -f "$W/prod.gro" ] && [ ! -f "$W/prod.xtc" ]; then
    echo "  $RUN: MDRUN GAGAL"; tail -6 "$W/prod.log" 2>/dev/null
    fail=$((fail+1)); consec=$((consec+1))
    # A transient fault burns the whole slice in seconds otherwise. On 2026-10-05 the GPU
    # briefly stopped being visible to mdrun -- "Compatible GPUs must have been found" --
    # and 23 runs failed back to back in minutes, each one instantly, leaving the worker
    # finished with a slice it had never actually run. Three in a row is not bad luck; it
    # is the machine, and the right response is to stop so the runs stay retryable rather
    # than racing through them.
    if [ "$consec" -ge 3 ]; then
      echo "  BERHENTI: $consec kegagalan mdrun berturut-turut -- kemungkinan GPU atau disk,"
      echo "            bukan ligan. Run yang belum dikerjakan tetap utuh untuk diulang."
      nvidia-smi --query-gpu=name,memory.used --format=csv,noheader 2>&1 | sed 's/^/            /'
      break
    fi
    continue
  fi
  consec=0
  el=$((SECONDS-t0))

  # descriptors, then the full trajectory can go
  if python "$ROOT/rung4_descriptors.py" --only "$RUN" --runs-dir "$ROOT/work" \
       --receptor "$ROOT/4ey7_rec_md.pdb" --out "$ROOT/out" $EXTRA_EXTRACT \
       >> "$ROOT/logs/extract.log" 2>&1; then
    if [ "$KEEP_REDUCED" = "1" ] && [ -f "$W/prod.xtc" ]; then
      # protein + ligand only, every 100 ps: ~50 MB instead of ~2 GB, so the option of
      # revisiting the trajectories later is not thrown away entirely
      # Build an explicit Protein-or-UNK union and select it BY NAME. Selecting group "0"
      # picks System, i.e. every atom -- and the file still shrinks, because -dt 100 cuts the
      # frame count, so the mistake looks like success. The union is always the last group
      # written, so read its name back rather than hard-coding a number that shifts with the
      # system's composition.
      printf '"Protein" | "UNK"\nq\n' | gmx make_ndx -f "$TPR" -o "$W/sel.ndx" > /dev/null 2>&1
      GRP=$(grep -o '^\[ [^]]* \]' "$W/sel.ndx" | tail -1 | tr -d '[] ')
      if [ -n "$GRP" ]; then
        gmx trjconv -s "$TPR" -f "$W/prod.xtc" -o "$ROOT/reduced/$RUN.xtc" \
            -dt 100 -n "$W/sel.ndx" -pbc mol <<< "$GRP" > /dev/null 2>&1 || true
        # A reduced trajectory without a matching topology is close to useless, and every
        # ligand has a different atom count, so each run needs its own reference structure.
        gmx trjconv -s "$TPR" -f "$W/prod.xtc" -o "$ROOT/reduced/$RUN.gro" \
            -dump 0 -n "$W/sel.ndx" -pbc mol <<< "$GRP" > /dev/null 2>&1 || true
      fi
    fi
    # The final checkpoint, always -- 2.3 MB against the 2 GB trajectory it replaces. It is
    # what keeps a shortened run reversible: with the tpr (kept on BRIN) a run can be picked
    # up later by
    #     gmx convert-tpr -s prod.tpr -extend <ps> -o next.tpr
    #     gmx mdrun -s next.tpr -cpi prod.cpt
    # so a decision to simulate 25 ns instead of 50 stops being a decision that cannot be
    # undone once the trajectory is gone. Without this, cutting the duration is a one-way bet.
    [ -f "$W/prod.cpt" ] && cp "$W/prod.cpt" "$ROOT/reduced/$RUN.cpt"
    rm -rf "$W"
    done_n=$((done_n+1))
    echo "  $RUN: selesai ${el}s, deskriptor OK, trajektori dihapus"
  else
    echo "  $RUN: EKSTRAKSI GAGAL -- trajektori DIPERTAHANKAN untuk diperiksa"
    tail -4 "$ROOT/logs/extract.log"; fail=$((fail+1))
  fi
done < "$SLICE"

echo
echo "selesai $done_n, dilewati $skip, gagal $fail"
echo "JSON di $ROOT/out: $(ls "$ROOT/out"/*.json 2>/dev/null | wc -l)"
echo "disk: $(df -h / | tail -1 | awk '{print $4" bebas"}')"
echo "berakhir : $(date)"
