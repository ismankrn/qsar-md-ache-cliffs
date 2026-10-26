#!/bin/bash
# Start exactly one worker on an instance, and refuse to start a second.
#
# This exists because of a real incident on 2026-10-23. Launching over ssh sometimes returns
# nothing and then drops the connection -- `setsid` detaches the worker, so the command
# SUCCEEDS precisely when the session that would have reported it dies. Reading that silence
# as failure and launching again put two workers on vast-3 and three on vast-4, all reading
# the same slice.txt, all picking the same run, and all running `mdrun -deffnm prod` inside
# one work directory. They overwrite each other's prod.xtc, prod.log and prod.cpt, and the
# trajectory is unusable. About 3.5 h of GPU was thrown away.
#
# So: never infer from a hung ssh that nothing started. Check with a FRESH connection, which
# is what this does before it launches anything.
#
# The pgrep pattern is bracketed. `pgrep -f vast_worker.sh` matches the ssh command line that
# carries the pattern itself, so an unbracketed check reports a worker on an idle machine --
# a false positive that has already wasted time twice in this project.
#
# Usage: bash md/launch_worker.sh <host> [host ...]
#        FORCE=1 bash md/launch_worker.sh <host>   # only after stopping the existing worker
SSHOPT=(-o BatchMode=yes -o ConnectTimeout=25 -o ServerAliveInterval=15 -o ServerAliveCountMax=3)
[ $# -eq 0 ] && { echo "pakai: bash md/launch_worker.sh vast-1 [vast-2 ...]"; exit 1; }

for H in "$@"; do
  echo "== $H"
  pre=$(ssh "${SSHOPT[@]}" "$H" 'echo "$(pgrep -cf "vast_[w]orker.sh") $(pgrep -cx gmx) $(ls -d /root/qsar/work/*/ 2>/dev/null | wc -l)"' 2>/dev/null | grep -vE "^(Welcome|Have fun|AI agents)" | tail -1)
  read -r W G D <<<"$pre"
  if [ -z "$W" ]; then echo "   tidak terjangkau -- tidak diluncurkan"; continue; fi
  if [ "$W" -gt 0 ] || [ "$G" -gt 0 ]; then
    if [ "${FORCE:-0}" != "1" ]; then
      echo "   SUDAH BERJALAN: worker=$W mdrun=$G work=$D -- TIDAK diluncurkan."
      echo "   Hentikan dulu kalau memang ingin memulai ulang (FORCE=1 hanya setelah itu)."
      continue
    fi
    echo "   worker=$W mdrun=$G sudah ada, FORCE=1 diberikan -- tetap diluncurkan"
  fi

  # The launch itself is expected to hang or drop; that is not an error. Verification comes
  # from a separate connection below, never from this command's output.
  ssh "${SSHOPT[@]}" "$H" 'cd /root/qsar && mkdir -p logs && \
     setsid nohup bash vast_worker.sh < /dev/null > logs/worker.$(date +%m%d-%H%M).log 2>&1 & \
     disown' > /dev/null 2>&1 &
  lpid=$!
  sleep 20
  kill "$lpid" 2>/dev/null

  post=$(ssh "${SSHOPT[@]}" "$H" 'cd /root/qsar; echo "$(pgrep -cf "vast_[w]orker.sh") $(pgrep -cx gmx) $(ls logs/worker.*.log 2>/dev/null | wc -l)"' 2>/dev/null | grep -vE "^(Welcome|Have fun|AI agents)" | tail -1)
  read -r W2 G2 L2 <<<"$post"
  if [ -z "$W2" ]; then echo "   peluncuran terkirim, verifikasi gagal terhubung -- CEK LAGI sebelum meluncurkan ulang"
  elif [ "$W2" -gt 0 ]; then echo "   berjalan: worker=$W2 mdrun=$G2 (log worker: $L2)"
  else echo "   TIDAK berjalan setelah peluncuran -- periksa logs/worker.*.log"; fi
done
