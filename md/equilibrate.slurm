#!/bin/bash
# Equilibration only (em -> nvt -> npt) for one batch. Production is a separate decision
# and a separate job; this work is identical whatever stage 11 concludes, so it can run now.
#
# `short`, one node, despite MaxNodes=1 there. medium-large would give 48 concurrent
# simulations against short's 12, but its queue put our jobs 47 h out while a new short
# job starts in 23 h: 6 h of work starting sooner beats 1.5 h of work starting a day
# later. Walltime does not affect the scheduler's estimate here -- it orders purely by
# priority, not backfill -- so asking for 20 h costs nothing.
#
# Production must still use medium-large: 50 ns at 16 cores takes 1.49 days, past short's
# 1-day cap, and medium-large's MaxNodes=4 is the only way to get real concurrency.
# Production is submitted only after equilibration finishes, so no two jobs ever touch the
# same run directory -- both phases use NBATCH=3, and one BATCH must have one owner.
# NOTE: these nodes are 2 sockets x 32 cores x 2 threads = 128 logical CPUs, but SLURM
# allocates physical cores: only 64 per node are schedulable. "8 tasks x 16 cpus = 128"
# is rejected with "Requested node configuration is not available". At 16 cores per
# simulation that is 4 simulations per node, not 8 -- which halves the concurrency the
# original budget assumed.
#SBATCH --job-name=md_equil
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --time=20:00:00
#SBATCH --output=/mgpfs/home/ikurniawan/qsar-md/logs/md_equil_b%a_%j.log

module load bioinformatics/gromacs/2023-plumed-mpi 2>/dev/null || \
module load bioinformatics/gromacs/2023.3-mpi
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

cd ~/qsar-md
echo "node  : $SLURM_JOB_NODELIST"
echo "batch : ${BATCH:-0} / ${NBATCH:-3}"
echo "mulai : $(date)"
echo

STOP_AFTER=npt BATCH=${BATCH:-0} NBATCH=${NBATCH:-3} \
  NTOMP=${SLURM_CPUS_PER_TASK:-16} MAXH=19 bash md/run_md.sh

echo
echo "selesai: $(date)"
