#!/usr/bin/env bash
# Progress report for a long-running stage.
#
#   ssh qsar 'cd qsar-md && bash scripts/progress.sh'       # stage 07
#   ssh qsar 'cd qsar-md && bash scripts/progress.sh 07'
#
# Safe to run at any time: it only reads files and never touches the job.
set -uo pipefail

STAGE="${1:-07}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY=/opt/miniforge3/envs/qsar-md/bin/python
[ -x "$PY" ] || PY=python3

case "$STAGE" in
  07)
    CSV="$ROOT/results/07_docking_scores.csv"
    LOG=/root/docking.log
    SESSION=dock
    TOTAL=1062
    LABEL="Stage 07 — production docking"
    ;;
  *)
    echo "unknown stage: $STAGE" >&2; exit 1 ;;
esac

echo "=============================================================="
echo " $LABEL"
echo "=============================================================="

# --- is it still running?
if tmux has-session -t "$SESSION" 2>/dev/null; then
    echo " status   : RUNNING   (tmux session '$SESSION')"
else
    echo " status   : NOT RUNNING   (tmux session '$SESSION' is gone)"
fi

if [ ! -f "$CSV" ]; then
    echo " no results yet — the first run has not finished"
    exit 0
fi

"$PY" - "$CSV" "$TOTAL" <<'PYEOF_INNER'
import sys, time, os
import pandas as pd
csv, total = sys.argv[1], int(sys.argv[2])
d = pd.read_csv(csv)
n = len(d)
age_h = (time.time() - os.path.getmtime(csv)) / 3600

med = d.seconds.median()
rate = 3600 / med if med else 0
remain = total - n
eta_h = remain * med / 3600 if med else 0

bar_w = 40
filled = int(bar_w * n / total)
bar = "#" * filled + "." * (bar_w - filled)

print(f" progress : [{bar}] {100*n/total:5.1f}%")
print(f" runs     : {n:,} of {total:,}   ({remain:,} left)")
print(f" speed    : {med:.0f} s per run  ->  {rate:.1f} runs/hour")
print(f" ETA      : {eta_h:.1f} hours  ({eta_h/24:.1f} days)")
print(f" last write: {age_h*60:.0f} minutes ago"
      + ("   <-- STALLED?" if age_h > 0.5 else ""))
print()
print(f" ligands done : {d.molecule_chembl_id.nunique()} unique")
print(f" best score   : {d.score_1.min():.2f} kcal/mol")
print(f" median score : {d.score_1.median():.2f} kcal/mol")
print(f" worst score  : {d.score_1.max():.2f} kcal/mol")

# seed agreement is the reason three seeds are run at all
full = d.groupby("molecule_chembl_id").filter(lambda g: len(g) >= 3)
if len(full):
    spread = full.groupby("molecule_chembl_id").score_1.agg(lambda s: s.max()-s.min())
    print()
    print(f" seed spread  : median {spread.median():.2f}, "
          f"max {spread.max():.2f} kcal/mol  ({len(spread)} ligands with 3 seeds)")
PYEOF_INNER

echo
FAILS=$(grep -c "FAIL" "$LOG" 2>/dev/null || echo 0)
echo " failures : $FAILS"
[ "$FAILS" -gt 0 ] && grep "FAIL" "$LOG" | tail -3 | sed 's/^/   /'

echo
echo " load     :$(uptime | sed 's/.*load average://')"
echo " disk     : $(df -h "$ROOT" | awk 'NR==2{print $4" free of "$2}')"
echo "=============================================================="
