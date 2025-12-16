#!/bin/bash
# submit_all.sh
# Submits one SBATCH job per (preset,cores) configuration using `run_single_template.slurm`.
# Adjust the PRESETS and CORES arrays below as desired.

set -euo pipefail

PRESETS=("small" "medium" "large")
CORES=(1 2 4 8)
CSV=tsp_preset_results.csv
RUNS=3

FIRST=1

for p in "${PRESETS[@]}"; do
  for c in "${CORES[@]}"; do
    if [ "$FIRST" -eq 1 ]; then
      OVERWRITE_FLAG="--overwrite"
      FIRST=0
    else
      OVERWRITE_FLAG=""
    fi

    echo "Submitting preset=$p cores=$c (overwrite=${OVERWRITE_FLAG})"
    # submit job requesting exactly c tasks
    sbatch --ntasks=${c} --export=ALL run_single_template.slurm ${p} ${c} ${CSV} ${RUNS} "${OVERWRITE_FLAG}"
    # small pause to avoid flooding the scheduler
    sleep 0.2
  done
done

echo "All jobs submitted. Monitor with 'squeue -u $(whoami)' or check output files."
