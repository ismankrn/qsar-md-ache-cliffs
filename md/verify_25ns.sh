#!/bin/bash
# Prove that what is running really is 25 ns, from the live process rather than from config.
#
# Reading the mdp would prove nothing here: prod.mdp still carries the campaign's original
# `nsteps = 25000000 ; 50 ns`, and so does every tpr built from it. What makes a run 25 ns is
# the `-nsteps 12500000` that vast_worker.sh puts on the mdrun command line, taken from the
# NSTEPS_PROD file that stage_vast.sh ships. If that file ever fails to arrive the worker
# drops the flag and runs the tpr's full 50 ns -- a trajectory that looks entirely healthy,
# takes twice as long, and is not comparable with the 217 development runs.
#
# So this reads three independent things: the file, the actual argv of the running mdrun, and
# the dt that mdrun itself reports. Finished runs are checked a fourth way, by frame count.
#
# Usage: bash md/verify_25ns.sh [host ...]
HOSTS=("$@"); [ ${#HOSTS[@]} -eq 0 ] && HOSTS=(vast-1 vast-2 vast-3 vast-4)
unproven=()
SSHOPT=(-o BatchMode=yes -o ConnectTimeout=20 -o ServerAliveInterval=15)
WANT=12500000
fail=0

for H in "${HOSTS[@]}"; do
  echo "== $H"
  out=$(ssh "${SSHOPT[@]}" "$H" '
    cd /root/qsar 2>/dev/null || { echo "NODIR"; exit; }
    echo "FILE $(cat NSTEPS_PROD 2>/dev/null || echo -)"
    for p in $(pgrep -x gmx); do
      tr "\0" " " < /proc/$p/cmdline | grep -o -- "-nsteps [0-9]*" | sed "s/^/ARGV /"
      R=$(basename "$(readlink /proc/$p/cwd)")
      echo "RUN $R"
      grep -E "^ *dt " "work/$R/prod.log" 2>/dev/null | head -1 | awk "{print \"DT \" \$3}"
    done
    for j in out/*.json; do
      [ -f "$j" ] || continue
      /venv/ana/bin/python -c "
import json,sys,os
d=json.load(open(sys.argv[1]))
print(\"JSON\", os.path.basename(sys.argv[1]), d[\"n_frames\"], d[\"ns_analysed\"])" "$j" 2>/dev/null
    done' 2>/dev/null | grep -vE "^(Welcome|Have fun|AI agents)")

  F=$(echo "$out" | sed -n 's/^FILE //p')
  if [ "$F" = "$WANT" ]; then echo "   NSTEPS_PROD = $F  OK"
  else echo "   NSTEPS_PROD = ${F:-tidak ada}  SALAH (harus $WANT)"; fail=1; fi

  A=$(echo "$out" | sed -n 's/^ARGV -nsteps //p' | head -1)
  D=$(echo "$out" | sed -n 's/^DT //p' | head -1)
  R=$(echo "$out" | sed -n 's/^RUN //p' | head -1)
  if [ -z "$A" ]; then
    # Not a pass. An instance still equilibrating has simply not been checked yet, and
    # reporting it as OK would be the exact failure this script exists to prevent.
    echo "   mdrun belum berjalan -- BELUM TERBUKTI, jalankan lagi setelah mdrun mulai"
    unproven+=("$H")
  elif [ "$A" = "$WANT" ]; then
    echo "   mdrun $R: -nsteps $A, dt $D ps -> $(awk -v a="$A" -v d="${D:-0.002}" 'BEGIN{printf "%g", a*d/1000}') ns  OK"
  else
    echo "   mdrun $R: -nsteps $A  SALAH"; fail=1
  fi

  echo "$out" | sed -n 's/^JSON //p' | while read -r name nf ns; do
    if [ "$nf" = "2001" ]; then echo "   $name: $nf frame, $ns ns teranalisis  OK"
    else echo "   $name: $nf frame, $ns ns  TIDAK SESUAI (harus 2001 / 20.01)"; fi
  done
  echo "$out" | sed -n 's/^JSON //p' | awk '$2!=2001{exit 1}' || fail=1
done

echo
if [ "$fail" -ne 0 ]; then
  echo "ADA YANG TIDAK SESUAI -- hentikan worker yang bersangkutan sebelum ia menulis JSON."; exit 1
elif [ ${#unproven[@]} -gt 0 ]; then
  echo "BELUM LENGKAP: ${#unproven[@]} instance belum menjalankan mdrun (${unproven[*]})."
  echo "Yang sudah berjalan terbukti 25 ns. Jalankan lagi untuk sisanya."; exit 2
else
  echo "SEMUA LOLOS: setiap simulasi berjalan 25 ns, dibuktikan dari argv mdrun yang hidup."
fi
