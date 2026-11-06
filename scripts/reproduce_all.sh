#!/usr/bin/env bash
# Run all four reproduction studies' canonical runners in sequence.
set -e
cd "$(dirname "$0")/.."
source .venv/bin/activate
for s in fig6-stochastic-trajectories fig4-nonresponder-dynamics fig2-ga-intervention fig3-generalizability; do
  echo "=================== $s ==================="
  python "workspace/studies/$s/sims/run.py"
done
echo "ALL REPRODUCTION STUDIES COMPLETE"
