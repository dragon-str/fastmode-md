#!/bin/bash
# Timed sweep for every fastmode system.
#
# The machine is shared, and fastmode only trusts a run when the 1-minute load
# average is below 2.0 at its start.  This script waits for that quiet window,
# then runs each system in turn.  Each audit writes its stdout to out/<name>.log
# so dashboard.py can show the in-flight candidate.
#
#   ./run_timed.sh          wait for load < 2.0, then sweep
#   ./run_timed.sh --now    start immediately (use only on a quiet machine)
set -u
cd "$(dirname "$0")"
# Set GMX and, if your force fields are not in the GROMACS installation, GMXLIB.
OUT=out
mkdir -p "$OUT"

quiet() { python3 -c 'import os,sys; sys.exit(0 if os.getloadavg()[0] < 2.0 else 1)'; }

if [ "${1:-}" != "--now" ]; then
    echo "waiting for load < 2.0 ..." | tee "$OUT/run_timed.status"
    until quiet; do sleep 60; done
fi
echo "machine quiet at $(date)" | tee -a "$OUT/run_timed.status"

run() {
    local name="$1"; shift
    echo "=== $name $(date) ===" | tee -a "$OUT/run_timed.status"
    python3 -u fastmode.py audit "$@" --strategy staged > "$OUT/$name.log" 2>&1
    echo "$name done $(date)" | tee -a "$OUT/run_timed.status"
}

run water_small  --gro systems/npt_3.0.gro --top systems/topol_3.0.top \
                 --out "$OUT/water_small"  --ref-ps 200
run water_medium --gro systems/npt_4.5.gro --top systems/topol_4.5.top \
                 --out "$OUT/water_medium" --ref-ps 200
run water_large  --gro systems/npt_6.0.gro --top systems/topol_6.0.top \
                 --out "$OUT/water_large"  --ref-ps 200
run villin       --gro systems/villin/eq.gro --top systems/villin/topol.top \
                 --out "$OUT/villin"       --ref-ps 200

echo "=== vsites $(date) ===" | tee -a "$OUT/run_timed.status"
python3 -u stage2.py sweep --ref-ps 200 > "$OUT/vsites.log" 2>&1
echo "ALL_DONE $(date)" | tee -a "$OUT/run_timed.status"
