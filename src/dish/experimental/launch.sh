#!/usr/bin/env bash
# Launch one experimental smoothing configuration and print its PID.
#
#   ./launch.sh                      # defaults: adapt_sigma=1, sigma=wall/4, CAD 12, p4
#   ./launch.sh 1 0.50 12 4          # sigma = wall/2
#   ./launch.sh 0 0.25 12 4          # uniform mollifier
#   ./launch.sh 1 0.25 16 4          # finer CAD quadrature for phi
#   ./launch.sh 1 0.25 12 4 50       # force rlam = 50
#
# args: adapt_sigma sigma_frac cad_order out_order [rlam_override]
set -euo pipefail

if (( $# > 5 )); then
    echo 'usage: launch.sh [adapt_sigma sigma_frac cad_order out_order rlam_override]' >&2
    exit 2
fi
AS=${1:-1}
SF=${2:-0.25}
CO=${3:-12}
OO=${4:-4}
RL=${5:-}

positive_number() {
    local number=$1
    [[ "$number" =~ ^[+]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$ ]] &&
        awk -v number="$number" 'BEGIN { exit !(number > 0 && number <= 1.7976931348623157e308) }'
}
if [[ ! "$AS" =~ ^[01]$ ]] || ! positive_number "$SF" ||
        [[ ! "$CO" =~ ^[1-9][0-9]*$ ]] || ! positive_number "$CO" ||
        [[ ! "$OO" =~ ^[1-9][0-9]*$ ]] || ! positive_number "$OO" ||
        { [[ -n "$RL" ]] && ! positive_number "$RL"; }; then
    echo 'Expected adapt_sigma 0 or 1, positive sigma/rlam, and positive integer orders.' >&2
    exit 2
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT_ROOT=$(cd -- "$SCRIPT_DIR/../../.." && pwd -P)
MESH_ROOT=${DISH_MESH_ROOT:-"$PROJECT_ROOT/data/dish"}
MESH_DIR="$MESH_ROOT/dish_fd023_m"
MATLAB_BIN=${MATLAB_BIN:-matlab}
if ! command -v "$MATLAB_BIN" >/dev/null 2>&1; then
    echo 'MATLAB was not found; set MATLAB_BIN to the MATLAB executable.' >&2
    exit 2
fi

TAG=$(awk -v as="$AS" -v sf="$SF" -v co="$CO" -v oo="$OO" \
    'BEGIN { printf "as%d_sig%03.0fw_cad%02d_p%02d", as, sf*100, co, oo }')
if [[ -n "$RL" ]]; then
    TAG="${TAG}_rlam$(awk -v rl="$RL" 'BEGIN { printf "%.17g", rl }')"
    CALL="run_smooth($AS,$SF,$CO,$OO,$RL)"
else
    CALL="run_smooth($AS,$SF,$CO,$OO)"
fi
LOG_DIR="$MESH_DIR/smooth/logs"
mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/$TAG.out"
# Pass the code path through the environment so quotes in paths remain safe.
export DISH_SMOOTH_CODE_DIR="$SCRIPT_DIR"
nohup "$MATLAB_BIN" -batch "addpath(getenv('DISH_SMOOTH_CODE_DIR')); $CALL" > "$LOG" 2>&1 &
PID=$!
echo "$PID" > "$LOG_DIR/$TAG.pid"
cat <<MSG
launched : $CALL
PID      : $PID
log      : $LOG
results  : $MESH_DIR/smooth/$TAG/
monitor  : tail -f "$LOG"
           ps -o pid=,etime=,pcpu=,rss= -p $PID
MSG
