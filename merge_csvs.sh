#!/usr/bin/env bash
# merge_csvs.sh
# Merge per-job CSV files (created by main when CSV_PER_JOB=1 or under SLURM)
# Usage: ./merge_csvs.sh combined.csv tsp_results.csv.job-*.csv

set -euo pipefail

# Behavior:
# - No args: auto-merge all files matching '*.job-*.csv' into
#   'tsp_sweep_results_combined.csv' (default output).
# - One arg: treat as output filename and auto-merge all '*.job-*.csv'.
# - Two or more args: first is output, remaining are input files (old behavior).

shopt -s nullglob

if [ "$#" -eq 0 ]; then
  out="tsp_sweep_results_combined.csv"
  # default to the exact prefix used by the run scripts
  inputs=( tsp_sweep_results.csv.job-*.csv )
elif [ "$#" -eq 1 ]; then
  out="$1"
  inputs=( tsp_sweep_results.csv.job-*.csv )
else
  out="$1"
  shift
  inputs=("$@")
fi

if [ ${#inputs[@]} -eq 0 ]; then
  # fallback: try find in current directory (helps when globbing behaved oddly)
  # fallback: try find for the exact prefix
  mapfile -t inputs < <(find . -maxdepth 1 -type f -name 'tsp_sweep_results.csv.job-*.csv' -printf '%f\n')
  if [ ${#inputs[@]} -eq 0 ]; then
    echo "No per-job CSV files found to merge (tried 'tsp_sweep_results.csv.job-*.csv' and 'find')." >&2
    echo "Run 'ls -1 *.job-*.csv' or 'ls -1 *job-*.csv' to see what's present." >&2
    exit 1
  fi
fi

echo "Found ${#inputs[@]} per-job CSV files to merge:" >&2
for f in "${inputs[@]}"; do echo "  $f" >&2; done

# write header from first readable file and append others (skip missing)
first=1
for f in "${inputs[@]}"; do
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

echo "Merged ${#inputs[@]} files into $out"
