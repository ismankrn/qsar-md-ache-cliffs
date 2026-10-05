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
FAILED=0
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
  # Never silence rsync. This is the only step that gets results off an instance that can be
  # reclaimed without warning, and it has already failed quietly once: three archives sat on
  # vast1 while the summary printed "0 baru" and looked like there was nothing to fetch. Vast
  # rate-limits rapid successive SSH connections, so a pull can fail for reasons that clear on
  # their own -- which is exactly why it must retry and then say so rather than shrug.
  pull_dir () {   # remote_subdir  local_dest  -> prints number of files transferred
    local sub=$1 dest=$2 try=0 rc=1 out=""
    while [ $try -lt 3 ]; do
      # --timeout and ConnectTimeout, because a hung transfer is worse than a failed one:
      # rsync with no I/O timeout blocks forever, and the monitor loop around it stops
      # dead. That is what silently froze the 10-minute poll for seventeen hours on
      # 2026-10-05 while its process stayed alive and looked healthy.
      out=$(rsync -a --update -i --timeout=120 \
            -e "ssh -o BatchMode=yes -o ConnectTimeout=25 -o ServerAliveInterval=15 -o ServerAliveCountMax=4" \
            "$H:/root/qsar/$sub/" "$dest/" 2>&1); rc=$?
      [ $rc -eq 0 ] && break
      try=$((try+1)); sleep 5
    done
    if [ $rc -ne 0 ]; then
      echo "  !! $H:$sub GAGAL setelah 3 percobaan (rc=$rc)" >&2
      echo "$out" | tail -3 | sed 's/^/     /' >&2
      FAILED=$((FAILED+1))
      echo 0; return
    fi
    echo "$out" | grep -c '^[<>]'
  }
  moved_j=$(pull_dir out "$DEST_JSON")
  moved_r=$(pull_dir reduced "$DEST_RED")
  after=$(ls "$DEST_JSON"/*.json 2>/dev/null | wc -l | tr -d ' ')
  fresh=$((after-before)); upd=$((moved_j-fresh)); [ $upd -lt 0 ] && upd=0
  slice=$(ssh -n -o BatchMode=yes "$H" 'wc -l < /root/qsar/slice.txt' 2>/dev/null)
  printf "  %-8s instance %s/%s selesai, +%s baru, %s diperbarui, %s berkas arsip\n" \
    "$H:" "$n_remote" "${slice:-?}" "$fresh" "$upd" "$moved_r"
  total_new=$((total_new + fresh))
done
echo "  ---"
# 217 phase-1 runs plus whatever replicates are in flight; counting the manifests beats a
# hard-coded 217 now that replicate JSONs land in the same directory.
WANT=$(( $(tail -n +2 "$ROOT/results/14_md_runs_phase1.csv" 2>/dev/null | grep -c .) \
       + $(tail -n +2 "$ROOT/results/14c_md_runs_repl_pilot.csv" 2>/dev/null | grep -c .) ))
echo "  JSON di laptop : $(ls "$DEST_JSON"/*.json 2>/dev/null | wc -l | tr -d ' ') / ${WANT:-217}"
echo "  arsip          : $(ls "$DEST_RED"/*.xtc 2>/dev/null | wc -l | tr -d ' ')  ($(du -sh "$DEST_RED" 2>/dev/null | cut -f1))"
echo "  baru kali ini  : $total_new"
# A pull that failed must not be readable as a pull that found nothing.
if [ "${FAILED:-0}" -gt 0 ]; then
  echo "  !! $FAILED transfer GAGAL -- hasil mungkin masih tertinggal di instance"
fi

# QC on every pull. A defective run used to be discoverable only by someone deciding to look,
# which is how two of them reached the repository; now nothing arrives without being checked.
# --new-only keeps this cheap: panels are re-rendered only where the archive is newer than the
# page, so a pull that brought nothing new costs a few seconds of reading JSONs.
#
# QC never fails the pull. Getting the data off a rented instance that can be reclaimed at any
# time is the job; reporting on it is secondary, and a missing matplotlib must not stop it.
# QC only when something actually arrived. A pull that finds nothing is two seconds of
# rsync; re-running QC over 80+ JSONs on top of that is what would make a 10-minute poll
# expensive. QC_ALWAYS=1 forces it.
if [ "${QC_ALWAYS:-0}" != "1" ] && [ "$total_new" -eq 0 ]; then
  exit 0
fi
PY=${QSAR_PYTHON:-/opt/anaconda3/envs/qsar-md/bin/python}
[ -x "$PY" ] || PY=$(command -v python3)
if [ -n "$PY" ] && "$PY" -c "import MDAnalysis, matplotlib" 2>/dev/null; then
  echo "  ---"
  "$PY" "$ROOT/scripts/rung4_qc_plots.py" --new-only 2>&1 | grep -v "^JSON records\|^panels already"
  rc=${PIPESTATUS[0]}
  [ "$rc" = "2" ] && echo "  >>> PERIKSA: ada run dengan cacat ekstraksi, lihat results/16c_qc_flags.csv"
else
  echo "  (QC dilewati: python dengan MDAnalysis+matplotlib tidak ditemukan)"
fi
