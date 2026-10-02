#!/usr/bin/env bash
#SBATCH --job-name=hirax-single-leaf
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --mail-user=yuguanw@princeton.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --chdir=/u/yw3384/hirax_maxwell/test
#SBATCH --output=/u/yw3384/hirax_maxwell/test/results/single_leaf_tip/polar-%j.out.log
#SBATCH --error=/u/yw3384/hirax_maxwell/test/results/single_leaf_tip/polar-%j.err.log

# Submit all default orders:
#   sbatch run_single_leaf_tip_polar.sh
# Submit selected orders:
#   sbatch run_single_leaf_tip_polar.sh "4 6 8"
set -euo pipefail

script_dir="/u/yw3384/hirax_maxwell/test"
matlab_bin="${MATLAB_BIN:-/usr/local/bin/matlab}"
orders="${1:-${HIRAX_SINGLE_LEAF_ORDERS:-}}"
cpus="${SLURM_CPUS_PER_TASK:?Missing Slurm CPU allocation}"

if [[ ! -r "$script_dir/run_single_leaf_tip.m" ]]; then
    printf 'Cannot read %s/run_single_leaf_tip.m.\n' "$script_dir" >&2
    exit 1
fi
if [[ ! -x "$matlab_bin" ]]; then
    printf 'MATLAB executable is unavailable: %s\n' "$matlab_bin" >&2
    exit 1
fi
if [[ -n "$orders" ]] && [[ ! "$orders" =~ ^[[:space:]]*[0-9]+([,[:space:]]+[0-9]+)*[[:space:]]*$ ]]; then
    printf 'Surface orders must be an increasing comma/space-separated integer list.\n' >&2
    exit 2
fi

if [[ -n "$orders" ]]; then
    export HIRAX_SINGLE_LEAF_ORDERS="$orders"
fi
export OMP_NUM_THREADS="$cpus"
export OMP_DYNAMIC=FALSE
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OPENBLAS_NUM_THREADS=1
export SRUN_CPUS_PER_TASK="$cpus"

mkdir -p -- "$script_dir/results/single_leaf_tip"
printf 'Started: %s\nJob: %s\nHost: %s\nCPUs: %s\nOrders: %s\n' \
    "$(date -Is)" "$SLURM_JOB_ID" "$(hostname -f)" "$cpus" \
    "${HIRAX_SINGLE_LEAF_ORDERS:-script default}"

matlab_command="maxNumCompThreads($cpus); set(groot,'defaultFigureVisible','off'); run('$script_dir/run_single_leaf_tip.m');"
srun --ntasks=1 --cpus-per-task="$cpus" \
    "$matlab_bin" -batch "$matlab_command"
