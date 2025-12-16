#!/bin/bash
# submit_sweep_from_template.sh
# Create and submit SBATCH jobs based on mpibatchfile.script, varying --ntasks and city count (--n=).
# Edits a temporary copy of the template for each (cores, n) pair and calls sbatch.

set -euo pipefail

TEMPLATE=mpibatchfile.script
CSV=tsp_sweep_results.csv
SEED=12345
RUNS=1
CORES=(1 2 4 8)
NS=(5 6 7 8 9)
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

for cores in "${CORES[@]}"; do
  for n in "${NS[@]}"; do
    TMP="./sbatch_tmp_${cores}t_${n}n.sh"
    cp "$TEMPLATE" "$TMP"

    # Ensure there is a SBATCH ntasks line; if present replace it, otherwise insert after first line
    if grep -q "^#SBATCH .*--ntasks" "$TMP"; then
      sed -i "s/^#SBATCH .*--ntasks=.*/#SBATCH --ntasks=${cores}/" "$TMP"
    else
      # insert after the first line
      sed -i "1a #SBATCH --ntasks=${cores}" "$TMP"
    fi

    # Replace any existing mpirun line with desired command
    # If cores==1 include --seq to get sequential timing
    OVERWRITE_FLAG=""
    if [ $FIRST -eq 1 ] && [ $OVERWRITE -eq 1 ]; then
      OVERWRITE_FLAG="--overwrite"
      FIRST=0
    fi

    if [ "$cores" -eq 1 ]; then
      CMD="mpirun ./main --n=${n} --runs=${RUNS} --csv=${CSV} ${OVERWRITE_FLAG} --seed=${SEED} --random --seq"
    else
      CMD="mpirun ./main --n=${n} --runs=${RUNS} --csv=${CSV} ${OVERWRITE_FLAG} --seed=${SEED} --random"
    fi

    # replace any line starting with mpirun or srun
    sed -i "/^\\s*\(mpirun\|srun\) /c\\${CMD}" "$TMP"

    echo "Prepared $TMP -> cores=$cores n=$n cmd: $CMD"
    if [ "$DRY_RUN" -eq 1 ]; then
      echo "DRY RUN: sbatch $TMP"
    else
      sbatch "$TMP"
    fi

    # small sleep to avoid submitting too quickly
    sleep 0.2
  done
done

echo "Done."
