#!/bin/bash
# Pull results, QC them, and bring a dead worker back.
#
# A Vast instance can reboot without warning: vast2 did on 2026-10-01 and its worker died
# with it, leaving the GPU idle and billing for 4.4 hours before anyone looked. The
# simulation state survived -- mdrun resumed from prod.cpt and was 27 ns in within five
# minutes -- so the only thing lost was the waiting. That is what this fixes: a worker that
# is dead while work remains gets started again, and the loss is bounded by one cycle
# instead of by how long it takes someone to notice.
#
# A restart is safe because the worker is idempotent: runs whose JSON exists are skipped,
# and an interrupted mdrun resumes from its checkpoint via -cpi.
#
# Usage: bash md/watch_vast.sh [interval_seconds] [host ...]
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
INTERVAL=${1:-10800}; shift 2>/dev/null
HOSTS=("$@"); [ ${#HOSTS[@]} -eq 0 ] && HOSTS=(vast1 vast2)
fails=0; restarts=0

probe () {   # host -> "RUNNING"|"STOPPED"|"UNREACHABLE"
  local out
  out=$(ssh -n -o BatchMode=yes -o ConnectTimeout=25 "$1" \
        'pgrep -f "vast_[w]orker" >/dev/null && echo RUNNING || echo STOPPED; echo RC_OK' 2>/dev/null)
  case "$out" in
    *RC_OK*RUNNING*|*RUNNING*RC_OK*) echo RUNNING ;;
    *RC_OK*)                         echo STOPPED ;;
    *)                               echo UNREACHABLE ;;
  esac
}

while :; do
  sleep "$INTERVAL"
  echo "[$(date '+%m-%d %H:%M')]"
  bash "$ROOT/md/pull_vast.sh" "${HOSTS[@]}" 2>&1
  alive=0; unreach=0
  for H in "${HOSTS[@]}"; do
    case "$(probe "$H")" in
      RUNNING) alive=$((alive+1)) ;;
      UNREACHABLE) unreach=$((unreach+1)); echo "   $H: tidak terjangkau" ;;
      STOPPED)
        # Count how many of THIS slice's runs have a result, not how many JSONs exist. out/
        # accumulates across slices -- vast2 held 45 after two slices and a replicate pilot
        # against a 28-line slice, so the old test read 45 >= 28 as "finished" and left a
        # dead worker alone for sixteen hours.
        n=$(ssh -n -o BatchMode=yes -o ConnectTimeout=25 "$H" \
            'c=0; while read -r r; do [ -n "$r" ] && [ -f "/root/qsar/out/$r.json" ] && c=$((c+1)); done < /root/qsar/slice.txt; echo $c' 2>/dev/null)
        s=$(ssh -n -o BatchMode=yes -o ConnectTimeout=25 "$H" 'grep -c . /root/qsar/slice.txt' 2>/dev/null)
        df=$(ssh -n -o BatchMode=yes -o ConnectTimeout=25 "$H" 'df -h / | tail -1 | awk "{print \$4}"' 2>/dev/null)
        # Never act on numbers that did not arrive: a failed ssh must not look like "no work done".
        if [ -z "$n" ] || [ -z "$s" ]; then
          echo "   $H: pekerja mati, tapi hitungan tidak terbaca -- tidak dijalankan ulang"
          unreach=$((unreach+1))
        elif [ "$n" -lt "$s" ]; then
          echo "   $H: pekerja MATI pada $n/$s selesai, disk $df -- dijalankan ulang"
          ssh -n -o BatchMode=yes -o ConnectTimeout=40 "$H" \
            "cd /root/qsar && mv -f logs/worker.log logs/worker.\$(date +%m%d-%H%M).log 2>/dev/null; NTOMP=\${NTOMP:-4} nohup bash vast_worker.sh > logs/worker.log 2>&1 &" 2>/dev/null
          sleep 20
          if [ "$(probe "$H")" = RUNNING ]; then
            restarts=$((restarts+1)); alive=$((alive+1))
            echo "   $H: hidup lagi (restart ke-$restarts)"
          else
            echo "   $H: GAGAL dijalankan ulang -- perlu diperiksa manual"
          fi
        else
          echo "   $H: slice tuntas ($n/$s)"
        fi
        ;;
    esac
  done
  echo "   hidup=$alive tak-terjangkau=$unreach restart-kumulatif=$restarts"
  if [ "$alive" -eq 0 ]; then
    if [ "$unreach" -gt 0 ]; then
      fails=$((fails+1))
      [ "$fails" -ge 8 ] && { echo "=== 8 siklus tanpa kontak, berhenti memantau ==="; break; }
    else
      echo "=== SEMUA SLICE TUNTAS ==="; break
    fi
  else fails=0; fi
done
bash "$ROOT/md/pull_vast.sh" "${HOSTS[@]}" 2>&1
