#!/bin/bash
# Install everything a rented instance needs to run one slice of the manifest: GROMACS with
# CUDA, and a python that can extract the rung-4 descriptors.
#
# A new Vast box is bare Ubuntu with the NVIDIA driver and nothing else. Two things about
# this that cost real time once each:
#
# 1. The GROMACS conda package brings no python. On 2026-10-23 all three instances ran
#    mdrun perfectly and then failed every extraction with "python: command not found",
#    which is survivable only because the worker keeps the trajectory when extraction
#    fails. Hence the `ana` env below. It is SEPARATE from `gmx` on purpose: a conda
#    transaction inside `gmx` can replace shared libraries that a running mdrun has mapped.
#
# 2. Env prefixes must be resolved, never assembled as $(conda info --base)/envs/<name>.
#    The Vast image's conda has base /root/miniforge3 but keeps envs in /venv, so every
#    hard-coded envs/ path silently misses -- which is why this script's own "gmx already
#    present" check never fired and it reinstalled GROMACS on an instance that had it.
#
# GROMACS 2023 and MDAnalysis 2.10.0 are pinned to match the 217 development runs. The CUDA
# build matters: the default conda-forge build is CPU-only and would run at a fraction of
# the speed while looking like it worked. A different MDAnalysis would put the test-split
# descriptors on a different code path, and their only purpose is comparability with the 217.
#
# Idempotent: safe to re-run on an instance that is mid-campaign.
#
# Usage: bash md/provision_vast.sh <host1> [host2 ...]
set -e
HOSTS=("$@")
[ ${#HOSTS[@]} -eq 0 ] && { echo "pakai: bash md/provision_vast.sh vast-1 [vast-2 ...]"; exit 1; }
SSHOPT=(-o BatchMode=yes -o ConnectTimeout=15 -o ServerAliveInterval=15
        -o ServerAliveCountMax=3)

for H in "${HOSTS[@]}"; do
  echo "== $H =="
  # No `ssh -n` here: it nulls stdin, which silently eats this heredoc and makes the whole
  # provisioning step report success having done nothing.
  ssh "${SSHOPT[@]}" "$H" 'bash -s' <<'REMOTE'
set -e
if [ ! -d /root/miniforge3 ]; then
  echo "   memasang miniforge..."
  curl -fsSL -o /tmp/mf.sh \
    https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash /tmp/mf.sh -b -p /root/miniforge3 > /tmp/mf.log 2>&1 \
    || { echo "   GAGAL memasang miniforge:"; tail -15 /tmp/mf.log; exit 1; }
  rm -f /tmp/mf.sh
fi
source /root/miniforge3/etc/profile.d/conda.sh
# Must never return non-zero. Under `set -e` an assignment takes the exit status of its
# command substitution, so `GMX=$(prefix_of gmx)` aborts the whole script the first time
# an env is absent -- which is exactly the case this function exists to detect. The
# failure is silent: provisioning returns 0 having printed nothing and installed nothing.
prefix_of() { ( conda activate "$1" 2>/dev/null && echo "$CONDA_PREFIX" ) || true; }

# --- GROMACS ---------------------------------------------------------------------------
GMX=$(prefix_of gmx)
if [ -n "$GMX" ] && [ -x "$GMX/bin/gmx" ]; then
  echo "   gmx sudah ada: $("$GMX/bin/gmx" --version 2>/dev/null | sed -n 's/^GROMACS version: *//p')"
else
  echo "   memasang GROMACS 2023 dengan CUDA..."
  # nompi_cuda: the single-node CUDA build. Without the selector conda-forge resolves the
  # CPU-only variant, which runs and is useless.
  conda create -y -n gmx -c conda-forge "gromacs=2023.*=nompi_cuda*" > /tmp/gmx.log 2>&1 \
    || conda create -y -n gmx -c conda-forge "gromacs=*=nompi_cuda*" > /tmp/gmx.log 2>&1 \
    || { echo "   GAGAL memasang GROMACS:"; tail -20 /tmp/gmx.log; exit 1; }
  GMX=$(prefix_of gmx)
fi
[ -n "$GMX" ] || { echo "   GAGAL: env gmx tidak dapat diaktifkan"; exit 1; }
V=$("$GMX/bin/gmx" --version 2>/dev/null | sed -n 's/^GROMACS version: *//p')
G=$("$GMX/bin/gmx" --version 2>/dev/null | grep -i 'GPU support' | head -1)
[ -n "$V" ] || { echo "   GAGAL: gmx tidak dapat dijalankan setelah pemasangan"; exit 1; }
echo "$G" | grep -qi CUDA || { echo "   GAGAL: build tanpa CUDA -- $G"; exit 1; }
echo "   gromacs  : $V   ($G)"

# --- interpreter for rung4_descriptors.py ----------------------------------------------
ANA=$(prefix_of ana)
if [ -z "$ANA" ] || ! "$ANA/bin/python" -c 'import MDAnalysis' 2>/dev/null; then
  echo "   memasang env analisis (MDAnalysis 2.10.0)..."
  conda create -y -n ana -c conda-forge python=3.11.16 \
      "mdanalysis=2.10.0" "numpy=1.26.4" "pandas=2.3.3" > /tmp/ana.log 2>&1 \
    || conda create -y -n ana -c conda-forge python=3.11 \
      "mdanalysis=2.10.0" > /tmp/ana.log 2>&1 \
    || { echo "   GAGAL memasang env analisis:"; tail -20 /tmp/ana.log; exit 1; }
  ANA=$(prefix_of ana)
fi
[ -n "$ANA" ] || { echo "   GAGAL: env ana tidak dapat diaktifkan"; exit 1; }
PV=$("$ANA/bin/python" -c 'import MDAnalysis,numpy,pandas as p;print(MDAnalysis.__version__,numpy.__version__,p.__version__)' 2>&1) \
  || { echo "   GAGAL: import tidak jalan -- $PV"; exit 1; }
echo "   analisis : MDAnalysis/numpy/pandas $PV"

# The worker runs `conda activate gmx`, so this directory is first on its PATH. Putting the
# shim here means a worker that is ALREADY RUNNING picks it up at its next extraction;
# teaching vast_worker.sh the path instead would require restarting the workers, and a
# restart treats the half-written prod.xtc of the run in flight as a finished trajectory,
# extracting a short run into a JSON that then permanently blocks the real one.
SHIM="$GMX/bin/python"
if [ -e "$SHIM" ] && ! grep -q 'envs-shim-qsar' "$SHIM" 2>/dev/null; then
  echo "   catatan: $SHIM sudah ada dan bukan shim kami -- dibiarkan"
else
  printf '#!/bin/bash\n# envs-shim-qsar\nexec %s/bin/python "$@"\n' "$ANA" > "$SHIM"
  chmod +x "$SHIM"
fi
PATH="$GMX/bin:$PATH"
python -c 'import MDAnalysis' || { echo "   GAGAL: shim python tidak berfungsi"; exit 1; }
echo "   python   : $(command -v python) -> ok"

# --- is this machine actually giving us the CPU it advertises? ------------------------
# vast-2 on 2026-10-23 advertised 32 CPUs and 98% reliability and delivered a single core's
# worth of throughput, 34x slower per core than its siblings. mdrun ran, the GPU sat at 0%
# waiting, and an energy minimisation that takes one minute elsewhere had not finished after
# twenty. Neither the reliability score nor the advertised core count showed it; only running
# something did. Two minutes here against twenty minutes of staging and an hour of confusion.
#
# The threshold is deliberately loose. Healthy instances measured 373, 390, 425 ms; the
# starved one 12,742 ms. Anything past 1500 ms is a different kind of machine, not a slower
# one, and nothing should be staged to it.
BENCH=$( { TIME_START=$(date +%s%N)
           "$ANA/bin/python" -c 's=0
for i in range(3000000): s+=i*i' >/dev/null 2>&1
           TIME_END=$(date +%s%N)
           echo $(( (TIME_END-TIME_START)/1000000 )); } )
echo "   cpu 1-core: ${BENCH} ms (sehat: ~400 ms)"
if [ "$BENCH" -gt 1500 ]; then
  echo "   DITOLAK: instance ini hanya mendapat sebagian kecil CPU yang diiklankan."
  echo "   Jangan distaging. Hapus dan sewa pengganti."
  exit 1
fi
REMOTE
done
echo
echo "selesai. Jalankan ulang worker di tiap instance."
