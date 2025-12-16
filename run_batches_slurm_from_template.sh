#!/bin/bash
#SBATCH -J TSP_BATCH
#SBATCH --ntasks=8
#SBATCH --export=ALL
#SBATCH --output=tsp-batch-%j.out
#SBATCH --time=0-02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE

# Derived from mpibatchfile.script: loads same MPI module then runs mpirun commands
module load openmpi/4.1.5/gcc/12.2.0

# Configuration (edit before submitting if desired)
BIN=main
CSV=tsp_preset_results.csv
RUNS=3
PRESETS=("small" "medium" "large")
# list of MPI ranks to test; set to the maximum you requested in --ntasks above
CORES=(1 2 4 8)

MPI_CMD="srun"

echo "Running TSP batch on Slurm node $(hostname)"

# Ensure binary exists (compile if needed)
if [ ! -x "$BIN" ]; then
  echo "Binary $BIN not found; compiling..."
  mpicc -O2 -o $BIN main.c -lm
fi

FIRST=1

for preset in "${PRESETS[@]}"; do
  echo "\n=== Preset: $preset ==="
  for cores in "${CORES[@]}"; do
    echo "--- Running preset=$preset on $cores MPI rank(s) (runs=$RUNS) ---"

    OVERWRITE_FLAG=""
    if [ "$FIRST" -eq 1 ]; then
      OVERWRITE_FLAG="--overwrite"
      FIRST=0
    fi

    if [ "$cores" -eq 1 ]; then
      echo "Running sequential baseline for preset=$preset"
      # srun will use one task within the allocation
      $MPI_CMD -n 1 ./$BIN --preset=$preset --runs=$RUNS --csv=$CSV $OVERWRITE_FLAG --seq
    else
      echo "Running parallel preset=$preset with $cores ranks"
      # request the desired number of tasks from Slurm's allocation
      $MPI_CMD -n $cores ./$BIN --preset=$preset --runs=$RUNS --csv=$CSV $OVERWRITE_FLAG
    fi

    echo "Completed preset=$preset cores=$cores"
    sleep 1
  done
done

echo "All preset runs completed. Results in: $CSV"
