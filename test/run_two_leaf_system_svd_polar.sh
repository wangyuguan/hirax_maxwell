#!/usr/bin/env bash
#SBATCH --job-name=hirax-two-leaf-svd
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --mail-user=yuguanw@princeton.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --chdir=/u/yw3384/hirax_maxwell/test
#SBATCH --output=/u/yw3384/hirax_maxwell/test/results/two_leaf_system_svd/polar-%j.out.log
#SBATCH --error=/u/yw3384/hirax_maxwell/test/results/two_leaf_system_svd/polar-%j.err.log

# Assemble the full two-leaf NRCCIE system in RAM and save its singular
# values (A itself is not stored):
#   sbatch run_two_leaf_system_svd_polar.sh
set -euo pipefail

script_dir="/u/yw3384/hirax_maxwell/test"
project_dir="$(dirname -- "$script_dir")"
matlab_bin="${MATLAB_BIN:-/usr/local/bin/matlab}"
cpus="${SLURM_CPUS_PER_TASK:?Missing Slurm CPU allocation}"

if [[ ! -x "$matlab_bin" ]]; then
    printf 'MATLAB executable is unavailable: %s\n' "$matlab_bin" >&2
    exit 1
fi
for required_file in \
    "$script_dir/run_two_leaf_system_svd.m" \
    "$project_dir/../fmm3dbie-hirax-dev/matlab/startup.m" \
    "$project_dir/../fmm3dbie-hirax-dev/matlab/fmm3dbie_routs.mexa64" \
    "$project_dir/../fmm3dbie-hirax-dev/FMM3D/matlab/fmm3d.mexa64" \
    "$project_dir/../chunkie/startup.m"; do
    if [[ ! -r "$required_file" ]]; then
        printf 'Required file is unavailable: %s\n' "$required_file" >&2
        exit 1
    fi
done

export NRCCIE_TWO_LEAF_ROWS_ONLY=0
export OMP_NUM_THREADS="$cpus"
export OMP_DYNAMIC=FALSE
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OPENBLAS_NUM_THREADS=1
export SRUN_CPUS_PER_TASK="$cpus"

mkdir -p -- "$script_dir/results/two_leaf_system_svd" "$project_dir/data"
printf 'Started: %s\nJob: %s\nHost: %s\nCPUs: %s\n' \
    "$(date -Is)" "$SLURM_JOB_ID" "$(hostname -f)" "$cpus"

matlab_command="maxNumCompThreads($cpus); set(groot,'defaultFigureVisible','off'); run('$script_dir/run_two_leaf_system_svd.m');"
srun --ntasks=1 --cpus-per-task="$cpus" \
    "$matlab_bin" -batch "$matlab_command"
