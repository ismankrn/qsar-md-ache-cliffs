"""Stage exactly the runs a host has been assigned, whatever its current state.

stage_vast.sh splits the manifest by stride and is the right tool at the start of a campaign.
It is the wrong tool after a rebalance: the stride decides which runs go where, and a
rebalance exists precisely because that assignment no longer matches the machines that are
alive. This takes the slice as given and makes the host match it.

Idempotent and additive. A system already on the instance is not re-sent, so re-staging a
host that has been working costs only the files it is actually missing.

Two things learned the hard way on 2026-10-23, both of which this encodes:

  * ONE ssh connection per host, not one per run. Thirteen back-to-back scp calls trip Vast's
    connection rate limiter and every one of them fails instantly -- and fails quietly, if
    stderr is redirected. A single `rsync --files-from` carries the whole set.

  * `-z`. system.gro and system.top are text and compress 6.7x. On the 0.14 MB/s uplink this
    turned 27 minutes of staging into 50 seconds.

Usage: python md/stage_runs.py <host> [host ...]      # reads /tmp/slice_<host>.txt
"""
import subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SSHOPT = ["-o","BatchMode=yes","-o","ConnectTimeout=25","-o","ServerAliveInterval=15"]
NOISE = ("Welcome to vast", "Have fun", "AI agents")

def clean(t): return "\n".join(l for l in t.splitlines() if not l.startswith(NOISE))
def ssh(h, cmd, timeout=180):
    r = subprocess.run(["ssh", *SSHOPT, h, cmd], capture_output=True, text=True, timeout=timeout)
    return clean(r.stdout)
def rsync(files_from, src, dest, timeout=1800):
    r = subprocess.run(["rsync","-az","--timeout=600",f"--files-from={files_from}",
                        "-e","ssh "+" ".join(SSHOPT), str(src), dest],
                       capture_output=True, text=True, timeout=timeout)
    return r.returncode, clean(r.stderr)

hosts = sys.argv[1:]
if not hosts: sys.exit("pakai: python md/stage_runs.py <host> [host ...]")

for h in hosts:
    slice_f = Path(f"/tmp/slice_{h}.txt")
    if not slice_f.exists(): print(f"{h}: {slice_f} tidak ada -- jalankan reslice.py dulu"); continue
    runs = [l.strip() for l in slice_f.read_text().split() if l.strip()]
    print(f"== {h}: {len(runs)} run")

    ssh(h, 'mkdir -p /root/qsar/{tpr,systems,mdp,work,out,reduced,logs}')

    # Infrastructure. NSTEPS_PROD rides along deliberately: without it the worker drops
    # -nsteps and runs the tpr's 25,000,000 steps, i.e. 50 ns instead of 25, which looks
    # healthy and is not comparable with the 217 development trajectories.
    infra = ["scripts/rung4_descriptors.py", "md/vast_worker.sh", "md/NSTEPS_PROD",
             "docking/4ey7_rec_md.pdb"]
    lf = Path(f"/tmp/infra_{h}.txt"); lf.write_text("\n".join(infra)+"\n")
    rc, err = rsync(lf, ROOT, f"{h}:/root/qsar/")
    if rc: print(f"   ! infra gagal: {err.splitlines()[-1] if err else rc}")
    ssh(h, 'cd /root/qsar && mv -f rung4_descriptors.py . 2>/dev/null; '
           'mv -f scripts/rung4_descriptors.py . 2>/dev/null; '
           'mv -f md/vast_worker.sh . 2>/dev/null; mv -f md/NSTEPS_PROD . 2>/dev/null; '
           'mv -f docking/4ey7_rec_md.pdb . 2>/dev/null; rmdir scripts md docking 2>/dev/null; true')

    lf = Path(f"/tmp/mdp_{h}.txt")
    lf.write_text("\n".join(f"mdp/{m}.mdp" for m in ("em","npt","prod"))+"\n")
    rc, err = rsync(lf, ROOT/"md", f"{h}:/root/qsar/")
    if rc: print(f"   ! mdp gagal: {err.splitlines()[-1] if err else rc}")

    have = set(ssh(h, 'ls -d /root/qsar/systems/*/ 2>/dev/null | xargs -r -n1 basename').split())
    need = [r for r in runs if r.rsplit("_r",1)[0] not in have]
    if need:
        lines = [f"systems/{r.rsplit('_r',1)[0]}/{f}" for r in need
                 for f in ("system.gro","system.top")]
        lf = Path(f"/tmp/sys_{h}.txt"); lf.write_text("\n".join(lines)+"\n")
        rc, err = rsync(lf, ROOT/"md", f"{h}:/root/qsar/")
        if rc: print(f"   ! sistem gagal: {err.splitlines()[-1] if err else rc}")
    print(f"   sistem: {len(have)} sudah ada, {len(need)} dikirim")

    # nvt.mdp carries the run's velocity seed, so it is per run and renamed on arrival.
    tmp = Path(f"/tmp/nvt_{h}"); subprocess.run(["rm","-rf",str(tmp)]); tmp.mkdir()
    for r in runs:
        src = ROOT/"md"/"runs"/r/"nvt.mdp"
        if src.exists(): (tmp/f"{r}.nvt.mdp").write_bytes(src.read_bytes())
    subprocess.run(["rsync","-az","--timeout=300","-e","ssh "+" ".join(SSHOPT),
                    f"{tmp}/", f"{h}:/root/qsar/mdp/"], capture_output=True, text=True)

    # Results already in hand go back, so the worker skips them instead of recomputing.
    back = [r for r in runs if (ROOT/"results/16_rung4_per_run"/f"{r}.json").exists()]
    if back:
        lf = Path(f"/tmp/json_{h}.txt")
        lf.write_text("\n".join(f"{r}.json" for r in back)+"\n")
        rsync(lf, ROOT/"results/16_rung4_per_run", f"{h}:/root/qsar/out/")

    subprocess.run(["scp", *SSHOPT, "-q", str(slice_f), f"{h}:/root/qsar/slice.txt"])

    chk = ssh(h, 'cd /root/qsar && echo "$(wc -l < slice.txt) $(ls -d systems/*/ 2>/dev/null | wc -l) '
                 '$(ls mdp/*.nvt.mdp 2>/dev/null | wc -l) $(cat NSTEPS_PROD 2>/dev/null) '
                 '$(ls out/*.json 2>/dev/null | wc -l) '
                 '$([ -f vast_worker.sh ] && [ -f rung4_descriptors.py ] && [ -f 4ey7_rec_md.pdb ] && echo y || echo n)"').split()
    if len(chk) == 6:
        sl, sy, nv, ns, js, ok = chk
        lig = len({r.rsplit('_r',1)[0] for r in runs})
        good = (sl == str(len(runs)) and int(sy) >= lig and int(nv) >= len(runs)
                and ns == "12500000" and ok == "y")
        print(f"   slice={sl} sistem={sy}(butuh>={lig}) nvt={nv} nsteps={ns} json={js} berkas={ok}"
              f"  -> {'SIAP' if good else 'BELUM SIAP'}")
    else:
        print(f"   verifikasi gagal: {chk}")
