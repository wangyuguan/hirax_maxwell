#!/usr/bin/env bash
# Run the single-leaf NRCCIE convergence test directly on Latte.
# Latte has no Slurm scheduler, so this script must be launched on Latte.
#
# Usage:
#   ./run_single_leaf_tip_latte.sh [CPU_COUNT] ["SURFACE_ORDERS"]
# Examples:
#   ./run_single_leaf_tip_latte.sh
#   ./run_single_leaf_tip_latte.sh 36 "4 6 8"
#   HIRAX_SINGLE_LEAF_GEOMETRY_ONLY=1 ./run_single_leaf_tip_latte.sh 4 "4"
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd -- "$script_dir/.." && pwd)"
matlab_bin="${MATLAB_BIN:-/usr/local/bin/matlab}"
workers="${1:-${HIRAX_CPUS:-36}}"
orders="${2:-${HIRAX_SINGLE_LEAF_ORDERS:-}}"

if [[ "$(hostname -s)" != latte* ]]; then
    printf 'This script is intended for Latte; current host is %s.\n' "$(hostname -f)" >&2
    exit 1
fi
if [[ ! "$workers" =~ ^[1-9][0-9]*$ ]] || (( workers > 72 )); then
    printf 'CPU_COUNT must be an integer from 1 to 72; received %s.\n' "$workers" >&2
    exit 2
fi
if [[ -n "$orders" ]] && [[ ! "$orders" =~ ^[[:space:]]*[0-9]+([,[:space:]]+[0-9]+)*[[:space:]]*$ ]]; then
    printf 'SURFACE_ORDERS must be an increasing comma/space-separated integer list.\n' >&2
    exit 2
fi
if [[ ! -x "$matlab_bin" ]]; then
    printf 'MATLAB executable is unavailable: %s\n' "$matlab_bin" >&2
    exit 1
fi
for required_file in \
    "$script_dir/run_single_leaf_tip.m" \
    "$project_dir/../fmm3dbie-hirax-dev/matlab/fmm3dbie_routs.mexa64" \
    "$project_dir/../fmm3dbie-hirax-dev/FMM3D/matlab/fmm3d.mexa64" \
    "$project_dir/../chunkie/startup.m"; do
    if [[ ! -r "$required_file" ]]; then
        printf 'Required file is unavailable: %s\n' "$required_file" >&2
        exit 1
    fi
done

if [[ -n "$orders" ]]; then
    export HIRAX_SINGLE_LEAF_ORDERS="$orders"
fi
export OMP_NUM_THREADS="$workers"
export OMP_DYNAMIC=FALSE
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OPENBLAS_NUM_THREADS=1

results_dir="$script_dir/results/single_leaf_tip"
mkdir -p -- "$results_dir"
run_tag="$(date +%Y%m%d_%H%M%S)_$$"
log_file="$results_dir/latte_${run_tag}.log"

printf 'Started: %s\nHost: %s\nCPUs: %s\nOrders: %s\nLog: %s\n' \
    "$(date -Is)" "$(hostname -f)" "$workers" \
    "${HIRAX_SINGLE_LEAF_ORDERS:-4 6 8 10 12}" "$log_file" | tee "$log_file"

matlab_command="maxNumCompThreads($workers); set(groot,'defaultFigureVisible','off'); run('$script_dir/run_single_leaf_tip.m');"
if command -v numactl >/dev/null 2>&1; then
    launcher=(numactl --interleave=all "$matlab_bin")
else
    launcher=("$matlab_bin")
fi

if "${launcher[@]}" -batch "$matlab_command" 2>&1 | tee -a "$log_file"; then
    status=0
else
    status=${PIPESTATUS[0]}
fi

printf '\nFinished: %s (exit %s)\nLog: %s\n' \
    "$(date -Is)" "$status" "$log_file" | tee -a "$log_file"
exit "$status"
