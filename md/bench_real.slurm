#!/bin/bash
# Benchmark on a real production system, and validate the equilibration chain while
# doing it. The 35.6 ns/day figure that the 5.5-day budget rests on was measured on a
# 92,935-atom system; the production systems are 96,260 atoms, 3.6% larger. Extrapolating
# would be guessing, and the batch layout depends on the answer.
#
# Equilibration runs first (em -> nvt -> npt) so the benchmark uses the actual prod.tpr,
# pressure coupling included -- not a stand-in built from unequilibrated coordinates.
#
# 16 cores, not 64. A 64-CPU request was given a start time 18 hours out behind 50 pending
# jobs in `short`; 16 cores backfills. And 16 is the number that matters -- it is the
# production configuration. The 32- and 64-core points were never decision-relevant.
#SBATCH --job-name=bench_real
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=/mgpfs/home/ikurniawan/qsar-md/logs/bench_real_%j.log

module load bioinformatics/gromacs/2023-plumed 2>/dev/null || \
module load bioinformatics/gromacs/2023.3-mpi
cd ~/qsar-md

RUN=${RUN:-$(tail -n +2 results/14_md_runs.csv | cut -d, -f1 | head -1)}
D="$HOME/qsar-md/md/runs/$RUN"
echo "run     : $RUN"
echo "node    : $(hostname), $SLURM_CPUS_PER_TASK core"
echo "gromacs : $(gmx --version 2>/dev/null | sed -n 's/^GROMACS version: *//p')"
echo "simd    : $(gmx --version 2>/dev/null | sed -n 's/^SIMD instructions: *//p')"
echo "atom    : $(head -2 $D/system.gro | tail -1 | tr -d ' ')"
echo "mulai   : $(date)"
echo

cd "$D"
stage () {   # name  coord  extra...
  local n=$1 c=$2; shift 2
  if [ -f "$n.gro" ]; then echo "  $n: sudah ada"; return 0; fi
  gmx grompp -f $n.mdp -c "$c" -p system.top -o $n.tpr -maxwarn 10 "$@" > grompp_$n.log 2>&1 \
    || { echo "  $n: GROMPP GAGAL"; grep -iA5 "ERROR\|Fatal" grompp_$n.log | head -12; return 1; }
  # no /usr/bin/time here: it is absent on the compute nodes and the wrapper took the
  # whole command down with it. SECONDS is a bash builtin and always available.
  local t0=$SECONDS
  gmx mdrun -deffnm $n -nt ${SLURM_CPUS_PER_TASK:-16} -pin on > mdrun_$n.log 2>&1 \
    || { echo "  $n: MDRUN GAGAL"; tail -15 mdrun_$n.log; tail -15 $n.log 2>/dev/null; return 1; }
  echo "  $n: OK  ($((SECONDS-t0)) s)"
}

echo "=== ekuilibrasi ==="
stage em  system.gro            || exit 1
stage nvt em.gro  -r em.gro     || exit 1
stage npt nvt.gro -r nvt.gro -t nvt.cpt || exit 1

echo
echo "=== suhu & densitas setelah NPT ==="
printf "Temperature\nDensity\nPotential\n0\n" | gmx energy -f npt.edr -o npt_terms.xvg 2>/dev/null \
  | sed -n "/Energy  *Average/,/^$/p"

echo
echo "=== throughput produksi pada sistem nyata ==="
gmx grompp -f prod.mdp -c npt.gro -p system.top -t npt.cpt -o prod.tpr -maxwarn 10 \
  > grompp_prod.log 2>&1 || { echo "GROMPP prod GAGAL"; exit 1; }
printf "  %-8s %12s %14s %18s\n" "core" "ns/day" "sim per 64c" "agregat ns/day"
printf "  %-8s %12s %14s %18s\n" "-------" "----------" "------------" "----------------"
for NT in 8 16; do
  rm -rf bench_$NT && mkdir bench_$NT && cd bench_$NT
  gmx mdrun -s ../prod.tpr -deffnm b -nsteps 40000 -resethway -noconfout \
      -nt $NT -pin on > b.out 2>&1
  P=$(awk '/^Performance:/{v=$2} END{printf "%.2f", v}' b.log 2>/dev/null)
  cd ..
  awk -v nt=$NT -v p="$P" 'BEGIN{n=int(64/nt); printf "  %-8s %12s %14d %18.2f\n", nt, p, n, p*n}'
done

echo
echo "selesai : $(date)"
