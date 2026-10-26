"""Redistribute the unfinished runs over the instances that actually exist now.

The stride in stage_vast.sh fixes a run to a slice index for the life of the campaign, which
is what makes the slices provably disjoint without a lock. The cost is that when an instance
dies its share does not move: after vast-1 was destroyed and vast-2 turned out to be CPU
starved, 25 of the 54 runs belonged to machines that were never going to run them, and no
surviving worker would ever pick them up.

Rebalancing is a stop-and-restart operation, never an edit of a live slice.txt. Appending to
a running worker's slice was tried on 2026-10-07 and read back correctly; truncating one on
2026-10-09 was not, and the loop went on to run two names from a buffer that no longer matched
the file, duplicating finished work. The loop must be stopped, the file replaced, the loop
started again -- a restart re-reads from the top and skips every run that already has a JSON.

One rule makes the restart nearly free: whatever a host has in flight stays on that host and
goes FIRST in its new slice. The worker then meets it immediately and `mdrun -cpi prod.cpt`
resumes it from its checkpoint instead of starting it over, so a rebalance costs minutes
rather than the hours already spent on the trajectory.

Writes nothing on its own: prints the plan and the slice files under /tmp for the driver to
push once the workers are down.

Usage: python md/reslice.py vast-3 vast-4 vast-5
"""
import json, subprocess, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT/"results/14m_md_runs_test.csv"
LOCAL_JSON = ROOT/"results/16_rung4_per_run"
SSHOPT = ["-o","BatchMode=yes","-o","ConnectTimeout=25","-o","ServerAliveInterval=15"]

hosts = sys.argv[1:]
if not hosts:
    sys.exit("pakai: python md/reslice.py <host> [host ...]")

def ssh(h, cmd):
    r = subprocess.run(["ssh", *SSHOPT, h, cmd], capture_output=True, text=True, timeout=120)
    return "\n".join(l for l in r.stdout.splitlines()
                     if not l.startswith(("Welcome to vast", "Have fun", "AI agents")))

runs = pd.read_csv(MANIFEST).run.tolist()

# Done = a JSON on the laptop. That is the only copy that survives an instance being
# destroyed, so it is the authority; an instance-side JSON that has not been pulled yet is
# counted by asking each host as well.
done = {r for r in runs if (LOCAL_JSON/f"{r}.json").exists()}
inflight = {}
for h in hosts:
    out = ssh(h, 'cd /root/qsar 2>/dev/null || exit; '
                 'ls out/*.json 2>/dev/null | xargs -r -n1 basename | sed "s/\\.json$//"; '
                 'echo "--"; '
                 'for p in $(pgrep -x gmx); do basename "$(readlink /proc/$p/cwd)"; done')
    part = out.split("--")
    done |= {x.strip() for x in part[0].split() if x.strip()}
    if len(part) > 1:
        f = [x.strip() for x in part[1].split() if x.strip()]
        if f: inflight[h] = f[0]

remaining = [r for r in runs if r not in done]
print(f"manifest   : {len(runs)} run")
print(f"selesai    : {len(done)}")
print(f"tersisa    : {len(remaining)}")
for h in hosts:
    print(f"  {h}: sedang mengerjakan {inflight.get(h, '(tidak ada)')}")

# Each host keeps its in-flight run, first in its slice, so the restart resumes from
# prod.cpt instead of repeating the trajectory.
slices = {h: ([inflight[h]] if inflight.get(h) in remaining else []) for h in hosts}
pool = [r for r in remaining if r not in {v for v in inflight.values()}]
for i, r in enumerate(pool):
    slices[hosts[i % len(hosts)]].append(r)

print()
for h in hosts:
    f = Path(f"/tmp/slice_{h}.txt")
    f.write_text("\n".join(slices[h]) + "\n")
    head = slices[h][0] if slices[h] else "-"
    tag = " (dilanjutkan dari checkpoint)" if inflight.get(h) == head else ""
    print(f"{h}: {len(slices[h])} run, mulai dari {head}{tag}  -> {f}")

tot = sum(len(v) for v in slices.values())
assert tot == len(remaining), f"BUG: {tot} != {len(remaining)}"
assert len({r for v in slices.values() for r in v}) == tot, "BUG: ada run yang muncul dua kali"
print(f"\n{tot} run terbagi, tidak ada yang ganda, tidak ada yang hilang.")
print(f"perkiraan kalender: {max(len(v) for v in slices.values()) * 2.11:.0f} jam")
