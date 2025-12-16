#!/bin/bash
#SBATCH --job-name=tsp-single
#SBATCH --output=tsp-single-%j.out
#SBATCH --error=tsp-single-%j.err
#SBATCH --time=0-02:00:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL

# Usage: submit with --ntasks set to desired number of MPI ranks:
# sbatch --ntasks=<cores> run_single_template.slurm <preset> <cores> [csv] [runs] [overwrite]

PRESET=${1:-small}
CORES_REQ=${2:-1}
CSV=${3:-tsp_preset_results.csv}
RUNS=${4:-3}
OVERWRITE_FLAG=${5:-}

echo "Starting job $SLURM_JOB_ID on $(hostname): preset=$PRESET cores=$CORES_REQ runs=$RUNS overwrite='$OVERWRITE_FLAG'"

module load openmpi/4.1.5/gcc/12.2.0

# compile if binary missing
if [ ! -x "main" ]; then
  echo "Binary 'main' not found; compiling..."
  mpicc -O2 -o main main.c -lm
fi

# Use srun to launch within the allocation. Ensure sbatch was invoked with --ntasks >= CORES_REQ.
srun -n ${CORES_REQ} ./main --preset=${PRESET} --runs=${RUNS} --csv=${CSV} ${OVERWRITE_FLAG}

echo "Job $SLURM_JOB_ID completed: preset=$PRESET cores=$CORES_REQ"
