# Traveling-Salesman-project
The goal is to parallelize some solutions to the traveling salesman problem

# usage instructions
In order to compile and run follow these steps:
1. Open terminal and navigate to the project directory
2. import the module for openmpi
   `module load openmpi/4.1.5/gcc/12.2.0`
3. Compile the code 
'''bash
    mpicc main.c -lm -o main
''' 
4. Run the script to submit jobs
   `./submit_sweep_from_template.sh --dry-run`  (to see what would be submitted without actually submitting)
   `./submit_sweep_from_template.sh` (to actually submit the jobs)

## Results will be appended to results.csv in the project directory