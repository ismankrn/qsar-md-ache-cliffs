#!/bin/bash
# One batch of MD runs. Submit three of these (BATCH=0,1,2 ... up to NBATCH-1) to fill
# the MaxJobs=3 allowance; resubmit the same BATCH to resume, since run_md.sh skips
# finished runs and restarts production from its checkpoint.
#
# 4 nodes x 64 schedulable cores = 256 cores; at 16 cores per simulation that is 16
# simulations per job, 48 across the three jobs MaxJobs allows. At the measured
# 33.61 ns/day, 379 trajectories of 50 ns come to about 12 days. That is the cluster
# ceiling, not a choice: medium-large is the only partition with MaxNodes > 1.
#
# 8 cores per simulation would give 96 concurrent and ~11.1 days, 6% better, but 50 ns
# then takes 2.78 days against a 3-day walltime -- five hours of margin. Not worth it,
# especially as both layouts need roughly four resubmissions either way.
#
# A job cannot finish its share in one go: 16 slots against ~126 runs per batch is eight
# chunks, and a 3-day walltime holds two. Resubmit the same BATCH to continue -- run_md.sh
# skips finished runs and restarts production from prod.cpt.
#
# medium-large, not short: 50 ns at ~35 ns/day is 1.4 days, past the 1-day cap on short.
# NOTE: these nodes are 2 sockets x 32 cores x 2 threads = 128 logical CPUs, but SLURM
# allocates physical cores: only 64 per node are schedulable. "8 tasks x 16 cpus = 128"
# is rejected with "Requested node configuration is not available". At 16 cores per
# simulation that is 4 simulations per node, not 8 -- which halves the concurrency the
# original budget assumed.
#SBATCH --job-name=md_prod
# medium-large only. It is the sole partition allowing more than one node per job
# (MaxNodes=4); short has 48 nodes but MaxNodes=1, and medium-small likewise. With
# MaxJobs=3 that puts the ceiling at 3 x 4 = 12 nodes, whatever else is free.
#SBATCH --partition=medium-large
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --time=3-00:00:00
#SBATCH --output=/mgpfs/home/ikurniawan/qsar-md/logs/md_prod_b%a_%j.log

module load bioinformatics/gromacs/2023-plumed-mpi 2>/dev/null || \
module load bioinformatics/gromacs/2023.3-mpi
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

cd ~/qsar-md
echo "node    : $SLURM_JOB_NODELIST"
echo "batch   : ${BATCH:-0} / ${NBATCH:-12}"
echo "gromacs : $(gmx_mpi --version 2>/dev/null | sed -n 's/^GROMACS version: *//p')"
echo "mulai   : $(date)"
echo

# -maxh is set just inside the 3-day walltime so mdrun checkpoints and exits cleanly
BATCH=${BATCH:-0} NBATCH=${NBATCH:-12} NTOMP=${SLURM_CPUS_PER_TASK:-16} MAXH=70 \
  bash md/run_md.sh

echo
echo "selesai : $(date)"
