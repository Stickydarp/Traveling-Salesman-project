#!/bin/bash
# submit_sweep_from_template.sh
# Create and submit SBATCH jobs based on mpibatchfile.script, varying --ntasks and city count (--n=).
# Edits a temporary copy of the template for each (cores, n) pair and calls sbatch.

set -euo pipefail

TEMPLATE=mpibatchfile.script
CSV=tsp_sweep_results.csv
# base seed; per-n seed = BASE_SEED + n
BASE_SEED=12400
RUNS=1
CORES=(1 2 4 8)
# include larger city counts up to 13 as requested
NS=(5 6 7 8 9 10 11)
DRY_RUN=0
OVERWRITE=1

usage() {
  cat <<EOF
Usage: $0 [--dry-run] [--csv name] [--seed s] [--runs r] [--cores "1 2 4"] [--ns "5 6 7"] [--no-overwrite]

Reads the SBATCH header and module lines from $TEMPLATE and creates a temporary SBATCH file per configuration.
It replaces the '#SBATCH --ntasks=...' line with the requested core count, and replaces the mpirun line
with a call to './main --n=<N> --runs=<RUNS> --csv=<CSV> [--seq]' (sequential baseline when cores=1 will include --seq).

Options:
  --dry-run         Print the sbatch commands instead of submitting
  --csv NAME        CSV filename (default: $CSV)
  --seed S          Seed for random coords (default: $SEED)
  --runs R          Number of runs per job (default: $RUNS)
  --cores "list"    Space-delimited core counts (default: ${CORES[*]})
  --ns "list"      Space-delimited city counts (default: ${NS[*]})
  --no-overwrite    Do not use --overwrite for the first job
EOF
}

# Parse args
while [ "$#" -gt 0 ]; do
  case "$1" in
    --dry-run) DRY_RUN=1; shift;;
    --csv) CSV="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --runs) RUNS="$2"; shift 2;;
    --cores) read -r -a CORES <<< "$2"; shift 2;;
    --ns) read -r -a NS <<< "$2"; shift 2;;
    --no-overwrite) OVERWRITE=0; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

if [ ! -f "$TEMPLATE" ]; then
  echo "Template $TEMPLATE not found in current directory." >&2
  exit 1
fi
FIRST=1

# Submit a sequential baseline for each n, then its parallel runs (cores > 1)
for n in "${NS[@]}"; do
  SEED=$((BASE_SEED + n))

  # submit sequential baseline (cores=1) first for this n
  cores=1
  TMP="./sbatch_tmp_${cores}t_${n}n.sh"
  cp "$TEMPLATE" "$TMP"
  if grep -q "^#SBATCH .*--ntasks" "$TMP"; then
    sed -i "s/^#SBATCH .*--ntasks=.*/#SBATCH --ntasks=${cores}/" "$TMP"
  else
    sed -i "1a #SBATCH --ntasks=${cores}" "$TMP"
  fi

  OVERWRITE_FLAG=""
  if [ $FIRST -eq 1 ] && [ $OVERWRITE -eq 1 ]; then
    OVERWRITE_FLAG="--overwrite"
    FIRST=0
  fi

  CMD_SEQ="mpirun ./main --n=${n} --runs=${RUNS} --csv=${CSV} ${OVERWRITE_FLAG} --seed=${SEED} --random --seq"
  sed -i "/^\\s*\(mpirun\|srun\) /c\\${CMD_SEQ}" "$TMP"
  echo "Prepared sequential $TMP -> cores=1 n=$n cmd: $CMD_SEQ"
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "DRY RUN: sbatch $TMP"
  else
    sbatch "$TMP"
  fi
  sleep 0.10

  # now submit parallel runs for this n
  for cores in "${CORES[@]}"; do
    if [ "$cores" -le 1 ]; then
      continue
    fi
    TMP="./sbatch_tmp_${cores}t_${n}n.sh"
    cp "$TEMPLATE" "$TMP"
    if grep -q "^#SBATCH .*--ntasks" "$TMP"; then
      sed -i "s/^#SBATCH .*--ntasks=.*/#SBATCH --ntasks=${cores}/" "$TMP"
    else
      sed -i "1a #SBATCH --ntasks=${cores}" "$TMP"
    fi

    CMD_PAR="mpirun ./main --n=${n} --runs=${RUNS} --csv=${CSV} ${OVERWRITE_FLAG} --seed=${SEED} --random"
    sed -i "/^\\s*\(mpirun\|srun\) /c\\${CMD_PAR}" "$TMP"

    echo "Prepared $TMP -> cores=$cores n=$n cmd: $CMD_PAR"
    if [ "$DRY_RUN" -eq 1 ]; then
      echo "DRY RUN: sbatch $TMP"
    else
      sbatch "$TMP"
    fi
    sleep 0.5
  done
done

echo "Done."
