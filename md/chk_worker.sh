#!/usr/bin/env bash
# How many vast_worker.sh instances are really running on this machine.
#
# `pgrep -fc vast_[w]orker` is the obvious test and it is wrong: the worker runs mdrun in a
# subshell, `( cd "$W" && gmx mdrun ... )`, and bash gives that subshell the PARENT's command
# line. So a single healthy worker shows up TWICE. On 2026-10-07 that nearly had me kill a
# perfectly good vast3, because a real duplicate -- two workers writing one run directory --
# had corrupted a run once before and 2 looked like the signature of it.
#
# The distinguishing feature is lineage: a top-level worker's parent is not itself a worker.
ps -eo pid,ppid,args | awk '
  $3 ~ /bash/ && $0 ~ /vast_worker\.sh/ { cmd[$1]=1; par[$1]=$2 }
  END { n=0; for (p in cmd) if (!(par[p] in cmd)) { n++; top=top " " p }
        printf "worker tingkat-atas: %d  (pid:%s)\n", n, top }'
# Real mdrun processes. Two traps here, both of which I walked into: `-e` overrides `-C`, so
# the gmx filter does nothing and every process on the box is listed; and a plain grep pattern
# matches the grep process's own command line. Hence the bracket, and no -e.
echo "mdrun            : $(ps -eo args= 2>/dev/null | grep -c 'gmx [m]drun')"
