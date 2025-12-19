#!/usr/bin/env bash
# merge_csvs.sh
# Merge per-job CSV files (created by main when CSV_PER_JOB=1 or under SLURM)
# Usage: ./merge_csvs.sh combined.csv tsp_results.csv.job-*.csv

set -euo pipefail
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <output.csv> <input1.csv> [input2.csv ...]"
  exit 2
fi
out="$1"
shift
# write header from first readable file and append others (skip missing)
first=1
for f in "$@"; do
  if [ ! -f "$f" ]; then
    echo "Warning: file '$f' not found, skipping" >&2
    continue
  fi
  if [ $first -eq 1 ]; then
    head -n 1 "$f" > "$out"
    tail -n +2 "$f" >> "$out"
    first=0
  else
    tail -n +2 "$f" >> "$out"
  fi
done

echo "Merged CSVs into $out"
