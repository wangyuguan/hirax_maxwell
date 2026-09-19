#!/bin/bash
#SBATCH --job-name=hirax-leaves-nrccie
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --time=48:00:00
#SBATCH --mail-user=yuguanw@princeton.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --chdir=/u/yw3384/hirax_maxwell
#SBATCH --output=/u/yw3384/hirax_maxwell/slurm_jobs/run_leaves_nrccie-%j.out
#SBATCH --error=/u/yw3384/hirax_maxwell/slurm_jobs/run_leaves_nrccie-%j.err

# Submit from hirax_maxwell: sbatch slurm_jobs/run_leaves_nrccie.sh
# Or from slurm_jobs:       sbatch run_leaves_nrccie.sh
# Supply --partition/--account to sbatch if required by your cluster.
# Optional: MATLAB_BIN=/path/to/matlab sbatch run_leaves_nrccie.sh
# Uses one MATLAB process with OpenMP threads matching the CPU allocation.
set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    printf 'Submit this script with sbatch.\n' >&2
    exit 1
fi

# Slurm sets the working directory using --chdir above.
if [[ ! -r run_leaves_nrccie.m ]]; then
    printf 'Cannot read run_leaves_nrccie.m in %s; check the Slurm --chdir setting.\n' "$PWD" >&2
    exit 1
fi

MATLAB_BIN="${MATLAB_BIN:-matlab}"
cpus="${SLURM_CPUS_PER_TASK:?Missing Slurm CPU allocation}"
export OMP_NUM_THREADS="$cpus"
export OMP_DYNAMIC=FALSE
export OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1
export SRUN_CPUS_PER_TASK="$cpus"

printf 'Job: %s | CPUs: %s | Directory: %s\n' "$SLURM_JOB_ID" "$cpus" "$PWD"
# Limit MATLAB numerical threads to the allocation and disable figure windows.
srun --ntasks=1 --cpus-per-task="$cpus" \
    "$MATLAB_BIN" -batch \
    "maxNumCompThreads(str2double(getenv('SLURM_CPUS_PER_TASK'))); set(groot,'defaultFigureVisible','off'); run_leaves_nrccie;"
