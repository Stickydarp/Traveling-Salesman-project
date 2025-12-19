#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>

typedef struct {
    int *path;       // contains the order of cities in this path ie path[0] = first city, path[1] = second city, ...
    double distance; // total distance of this path is zero in not evaluated yet
} coordinateSet;

typedef struct {
    int *xcoords; 
    int *ycoords;
    int numCities;
    int *bestPath;      
    double bestDistance;
    double *costMatrix; /* flattened n*n matrix: costMatrix[i*n + j] */
} travelingSalesman;

/* initialize the struct (does not copy arrays; just stores pointers) */
void travelingSalesman_init(travelingSalesman *ts, int *x, int *y, int n) {
    ts->xcoords = x;
    ts->ycoords = y;
    ts->numCities = n;
    ts->bestPath = NULL;
    ts->bestDistance = DBL_MAX;
    ts->costMatrix = NULL;
}

/* build an n*n flattened cost matrix (Euclidean distances) */
void travelingSalesman_buildCostMatrix(travelingSalesman *ts) {
    int n = ts->numCities;
    ts->costMatrix = malloc(sizeof(double) * n * n);
    if (!ts->costMatrix) {
        fprintf(stderr, "Failed to allocate cost matrix\n");
        return;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double dx = ts->xcoords[i] - ts->xcoords[j];
            double dy = ts->ycoords[i] - ts->ycoords[j];
            ts->costMatrix[i * n + j] = sqrt(dx * dx + dy * dy);
        }
    }
}

void travelingSalesman_free(travelingSalesman *ts) {
    if (ts->bestPath) {
        free(ts->bestPath);
        ts->bestPath = NULL;
    }
    if (ts->costMatrix) {
        free(ts->costMatrix);
        ts->costMatrix = NULL;
    }
}

// evaluate the total distance of a given path using costMatrix if available
double evaluatePath(travelingSalesman *ts, int *path) {
    double totalDistance = 0.0;
    int n = ts->numCities;
    if (ts->costMatrix) {
        for (int i = 0; i < n; i++) {
            int a = path[i];
            int b = path[(i + 1) % n];
            totalDistance += ts->costMatrix[a * n + b];
        }
        return totalDistance;
    }
    /* fallback to on-the-fly calculation */
    for (int i = 0; i < n; i++) {
        int cityA = path[i];
        int cityB = path[(i + 1) % n];
        double dx = ts->xcoords[cityA] - ts->xcoords[cityB];
        double dy = ts->ycoords[cityA] - ts->ycoords[cityB];
        totalDistance += sqrt(dx * dx + dy * dy);
    }
    return totalDistance;
}


/* prototypes for helper functions used by RDFS */
int *remove_index_copy(const int *arr, int len, int idx);
void remove_index_inplace(int *arr, int *len, int idx);

/* calls itself recursively: add remaining_len and path_len to signature */
coordinateSet RDFS(travelingSalesman *ts, int *remaining, int remaining_len, int *path, int path_len) {
    double minDist = DBL_MAX;
    coordinateSet bestSet;
    bestSet.path = NULL;
    bestSet.distance = DBL_MAX;
    //traversing the tree 
    for (int remIndex = 0; remIndex < remaining_len; remIndex++) {
        int city = remaining[remIndex];
        int *new_remaining = remove_index_copy(remaining, remaining_len, remIndex);

        int *new_path = malloc(sizeof(int) * (path_len + 1));
        memcpy(new_path, path, sizeof(int) * path_len);
        new_path[path_len] = city;

        coordinateSet currentSet = RDFS(ts, new_remaining, remaining_len - 1, new_path, path_len + 1);

        /* new_remaining used only by the recursive call, free it now */
        free(new_remaining);

        /* DO NOT free new_path here — currentSet.path may point to it.
           free currentSet.path only when the candidate is not chosen. */
        
        if (currentSet.distance < minDist) {
            if (bestSet.path) {    /* only free if previously allocated */
                free(bestSet.path);
            }
            bestSet = currentSet;
            minDist = currentSet.distance;
        } else {
            /* avoid leaking currentSet.path when not chosen */
            if (currentSet.path) free(currentSet.path);
        }
    }
    if(remaining_len==0){
        bestSet.distance=evaluatePath(ts,path);
        bestSet.path=path;
    }
    return bestSet;
}


/* 
uses DFS to find the shortest path in sequential brute force manner
*/
coordinateSet travelingSalesman_SEQ_bruteForce(travelingSalesman *ts) {
    int n = ts->numCities;
    coordinateSet bestPath;
    bestPath.path = NULL;
    bestPath.distance = DBL_MAX;

    if (n <= 0) return bestPath;

    /* build cost matrix for fast distance lookups */
    travelingSalesman_buildCostMatrix(ts);

    if (!ts->costMatrix) {
        fprintf(stderr, "Failed to build cost matrix\n");
        return bestPath;
    }

    /* prepare remaining cities (exclude start = 0) */
    int remaining_len = (n > 1) ? n - 1 : 0;
    int *remaining = NULL;
    if (remaining_len > 0) {
        remaining = malloc(sizeof(int) * remaining_len);
        if (!remaining) {
            fprintf(stderr, "Allocation failed\n");
            return bestPath;
        }
        for (int i = 0; i < remaining_len; ++i) remaining[i] = i + 1;
    }

    /* initial path starts at city 0 */
    int *initial_path = malloc(sizeof(int) * 1);
    if (!initial_path) {
        free(remaining);
        return bestPath;
    }
    initial_path[0] = 0;

    /* run recursive DFS */
    bestPath = RDFS(ts, remaining ? remaining : NULL, remaining_len, initial_path, 1);

    /* clean up: free remaining and initial_path (if not returned) */
    free(remaining);
    if (!(bestPath.path == initial_path)) {
        free(initial_path);
    }
    return bestPath;
}


//parallel brute force version
//each process will explore a subset of the permutations
coordinateSet travelingSalesman_PAR_bruteForce(travelingSalesman *ts, int mpi_rank, int mpi_size) {
    int n = ts->numCities;
    coordinateSet bestPath;
    bestPath.path = NULL;
    bestPath.distance = DBL_MAX;

    coordinateSet localBest;
    localBest.path = NULL;
    localBest.distance = DBL_MAX;

    if (n <= 0) return bestPath;

    /* build cost matrix on rank 0 and broadcast it to all ranks */
    if (mpi_rank == 0) {
        travelingSalesman_buildCostMatrix(ts);
        if (!ts->costMatrix) {
            fprintf(stderr, "Failed to build cost matrix on root\n");
            return bestPath;
        }
        MPI_Bcast(ts->costMatrix, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        ts->costMatrix = malloc(sizeof(double) * n * n);
        if (!ts->costMatrix) {
            fprintf(stderr, "Failed to allocate cost matrix on rank %d\n", mpi_rank);
            return bestPath;
        }
        MPI_Bcast(ts->costMatrix, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    /* Distribute second-city choices across ranks.
       Each branch fixes city `c` as the second city (after 0) and explores the subtree. */
    for (int city = 1; city < n; ++city) {
        if (((city - 1) % mpi_size) != mpi_rank) continue; //skips the loop if not this rank's turn

        /* build remaining array: all cities except 0 and `city` */
        int remaining_len = (n > 2) ? (n - 2) : 0;
        int *remaining = NULL;
        if (remaining_len > 0) {
            remaining = malloc(sizeof(int) * remaining_len);
            if (!remaining) {
                fprintf(stderr, "Allocation failed on rank %d\n", mpi_rank);
                /* free localBest.path if allocated and return */
                if (localBest.path) free(localBest.path);
                return bestPath;
            }
            int idx = 0;
            for (int j = 1; j < n; ++j) {
                if (j == city) continue;
                remaining[idx++] = j;
            }
        }

        /* initial path = [0, city] */
        int *initial_path = malloc(sizeof(int) * 2);
        if (!initial_path) {
            free(remaining);
            if (localBest.path) free(localBest.path);
            return bestPath;
        }
        initial_path[0] = 0;
        initial_path[1] = city;

        /* explore subtree starting from this fixed second city */
        coordinateSet branchBest = RDFS(ts, remaining, remaining_len, initial_path, 2);

        /* branchBest.path may point to initial_path (or to deeper allocated arrays).
           If branchBest is better than localBest, adopt it; otherwise free branchBest.path. */
        if (branchBest.distance < localBest.distance) {
            if (localBest.path) free(localBest.path);
            localBest = branchBest;
        } else {
            if (branchBest.path) free(branchBest.path);
        }

        free(remaining);
        /* do not free initial_path here if branchBest.path points to it (it was moved) */
    }

    /* Now gather local best distances on root to determine global best rank */
    double localDist = localBest.distance;
    double *allDists = NULL;
    if (mpi_rank == 0) {
        allDists = malloc(sizeof(double) * mpi_size);
        if (!allDists) {
            fprintf(stderr, "Root failed to allocate allDists\n");
            if (localBest.path) free(localBest.path);
            return bestPath;
        }
    }
    MPI_Gather(&localDist, 1, MPI_DOUBLE, allDists, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int bestRank = -1;
    double globalBestDist = DBL_MAX;
    if (mpi_rank == 0) {
        /* find rank with minimum distance */
        for (int r = 0; r < mpi_size; ++r) {
            if (allDists[r] < globalBestDist) {
                globalBestDist = allDists[r];
                bestRank = r;
            }
        }
    }

    /* gather full paths (each path is n ints). processes without a path send -1-filled arrays */
    int *localPathBuf = malloc(sizeof(int) * n);
    if (!localPathBuf) {
        fprintf(stderr, "Rank %d failed to allocate localPathBuf\n", mpi_rank);
        if (localBest.path) free(localBest.path);
        if (allDists) free(allDists);
        return bestPath;
    }
    if (localBest.path) {
        for (int i = 0; i < n; ++i) localPathBuf[i] = localBest.path[i];
    } else {
        for (int i = 0; i < n; ++i) localPathBuf[i] = -1;
    }

    int *allPaths = NULL;
    if (mpi_rank == 0) {
        allPaths = malloc(sizeof(int) * n * mpi_size);
        if (!allPaths) {
            fprintf(stderr, "Root failed to allocate allPaths\n");
            free(localPathBuf);
            if (localBest.path) free(localBest.path);
            if (allDists) free(allDists);
            return bestPath;
        }
    }

    MPI_Gather(localPathBuf, n, MPI_INT, allPaths, n, MPI_INT, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        if (bestRank >= 0 && allPaths) {
            /* copy best path into return structure */
            bestPath.path = malloc(sizeof(int) * n);
            if (!bestPath.path) {
                fprintf(stderr, "Root failed to allocate bestPath.path\n");
            } else {
                for (int i = 0; i < n; ++i) {
                    bestPath.path[i] = allPaths[bestRank * n + i];
                }
                bestPath.distance = globalBestDist;
            }
        }
    }

    /* cleanup local allocations */
    free(localPathBuf);
    if (allPaths) free(allPaths);
    if (allDists) free(allDists);
    if (localBest.path) free(localBest.path);

    /* Note: non-root ranks return bestPath.path == NULL, distance DBL_MAX.
       Root returns the global best. */
    return bestPath;
}

//creates a copy of the given list while removing a specific index
int *remove_index_copy(const int *arr, int len, int idx) {
    if (idx < 0 || idx >= len) return NULL;
    int *out = malloc(sizeof(int) * (len - 1));
    if (!out) return NULL;
    for (int i = 0, j = 0; i < len; ++i) {
        if (i == idx) continue;
        out[j++] = arr[i];
    }
    return out;
}

void remove_index_inplace(int *arr, int *len, int idx) {
    if (idx < 0 || idx >= *len) return;
    memmove(&arr[idx], &arr[idx + 1], (*len - idx - 1) * sizeof(int));
    (*len)--;
}

int main(int argc, char *argv[]) {
    int rank, size;
    int force_seq = 0;
    int use_random = 0;
    int preset = 0; /* 0=small(default), 1=medium, 2=large */
    int custom_n = 0;
    unsigned int rand_seed = (unsigned int)time(NULL);

    /* new batch/CSV options */
    char csv_filename[512] = "tsp_results.csv";
    int runs_per_n = 1;
    int min_n = 0, max_n = 0, step_n = 1;
    int overwrite_csv = 0;

    /* parse simple command-line flags before MPI_Init */
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--seq") == 0) force_seq = 1;
        else if (strncmp(argv[i], "--n=", 4) == 0) custom_n = atoi(argv[i] + 4);
        else if (strcmp(argv[i], "--random") == 0) use_random = 1;
        else if (strncmp(argv[i], "--seed=", 7) == 0) rand_seed = (unsigned int)atoi(argv[i] + 7);
        else if (strcmp(argv[i], "--preset=medium") == 0) preset = 1;
        else if (strcmp(argv[i], "--preset=large") == 0) preset = 2;
        else if (strncmp(argv[i], "--csv=", 6) == 0) strncpy(csv_filename, argv[i] + 6, sizeof(csv_filename) - 1);
        else if (strncmp(argv[i], "--runs=", 7) == 0) runs_per_n = atoi(argv[i] + 7);
        else if (strncmp(argv[i], "--min-n=", 8) == 0) min_n = atoi(argv[i] + 8);
        else if (strncmp(argv[i], "--max-n=", 8) == 0) max_n = atoi(argv[i] + 8);
        else if (strncmp(argv[i], "--step=", 7) == 0) step_n = atoi(argv[i] + 7);
        else if (strcmp(argv[i], "--overwrite") == 0) overwrite_csv = 1;
    }

    MPI_Init(&argc, &argv);                 /* Initialize the MPI environment */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   /* Get the rank of the process */
    MPI_Comm_size(MPI_COMM_WORLD, &size);   /* Get the number of processes */

    if (rank == 0) printf("Comparing sequential and parallel TSP on %d MPI ranks\n", size);

    /* default small dataset */
    int small_n = 5;
    int small_xs_stack[5] = {0,1,2,3,5};
    int small_ys_stack[5] = {0,1,2,3,8};

    /* medium and large presets (increase as desired; beware factorial growth) */
    int medium_n = 8;
    int medium_xs_stack[8] = {0,2,5,9,4,7,3,8};
    int medium_ys_stack[8] = {0,1,6,2,9,5,8,3};

    int large_n = 10;
    int large_xs_stack[10] = {0,10,20,30,15,25,5,35,12,22};
    int large_ys_stack[10] = {0,5,15,8,20,2,25,10,18,12};

    /* decide iteration range */
    if (min_n > 0 && max_n > 0 && max_n < min_n) {
        if (rank == 0) fprintf(stderr, "Invalid min/max n range\n");
        MPI_Finalize();
        return 1;
    }

    if (custom_n > 0) {
        if (min_n == 0 && max_n == 0) {
            min_n = max_n = custom_n;
        }
        use_random = 1; /* custom n implies generating coords unless user supplies presets */
    }

    if (min_n == 0 && max_n == 0) {
        /* No range provided, use preset or single n from preset */
        if (preset == 1) { min_n = max_n = medium_n; }
        else if (preset == 2) { min_n = max_n = large_n; }
        else if (custom_n > 0) { min_n = max_n = custom_n; }
        else { min_n = max_n = small_n; }
    }

    if (step_n <= 0) step_n = 1;

    /* CSV file handling on rank 0
       To avoid multiple parallel jobs concurrently appendig to the file*/
    FILE *csvf = NULL;
    char final_csv[1024];
    if (rank == 0) {
        const char *per_job = getenv("CSV_PER_JOB");
        const char *slurm_job = getenv("SLURM_JOB_ID");
        if (per_job != NULL || slurm_job != NULL) {
            /* build unique filename: base + .job-<jobid or host>-<pid>.csv */
            const char *idpart = slurm_job ? slurm_job : "local";
            char hostbuf[128] = {0};
            gethostname(hostbuf, sizeof(hostbuf));
            snprintf(final_csv, sizeof(final_csv), "%s.job-%s-%s-%d.csv", csv_filename, idpart, hostbuf, (int)getpid());
            if (overwrite_csv) csvf = fopen(final_csv, "w");
            else csvf = fopen(final_csv, "a");
        } else {
            /* default behavior: append to provided CSV filename */
            strncpy(final_csv, csv_filename, sizeof(final_csv)-1);
            final_csv[sizeof(final_csv)-1] = '\0';
            if (overwrite_csv) {
                csvf = fopen(final_csv, "w");
            } else {
                FILE *tf = fopen(final_csv, "r");
                if (tf) { fclose(tf); csvf = fopen(final_csv, "a"); }
                else csvf = fopen(final_csv, "w");
            }
        }
        if (!csvf) {
            fprintf(stderr, "Failed to open CSV file '%s' for writing\n", final_csv);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        if (overwrite_csv || ftell(csvf) == 0) {
            fprintf(csvf, "n,mpi_ranks,run,seed,seq_time_s,par_time_s,best_distance\n");
            fflush(csvf);
        }
        if (per_job != NULL || slurm_job != NULL) {
            printf("Writing per-job CSV output to '%s'\n", final_csv);
        } else {
            printf("Writing CSV output to '%s'\n", final_csv);
        }
    }

    /* iterate over requested n values and perform runs_per_n runs */
    for (int n = min_n; n <= max_n; n += step_n) {
        if (rank == 0) {
            printf("Starting tests for n=%d (runs=%d)\n", n, runs_per_n);
            if (n <= 12) {
                printf("Note: brute-force TSP scales as n! — n=%d may be expensive\n", n);
            }
        }

        for (int run = 1; run <= runs_per_n; ++run) {
            int *xs = NULL;
            int *ys = NULL;
            int allocated_coords = 0;
            unsigned int seed = rand_seed + run + n * 13;

            if (!use_random) {
                /* use preset arrays when sizes match; otherwise fall back to random */
                if (n == small_n) {
                    xs = small_xs_stack; ys = small_ys_stack;
                } else if (n == medium_n) {
                    xs = medium_xs_stack; ys = medium_ys_stack;
                } else if (n == large_n) {
                    xs = large_xs_stack; ys = large_ys_stack;
                } else {
                    /* fall back to random when preset not available */
                    allocated_coords = 1;
                    xs = malloc(sizeof(int) * n);
                    ys = malloc(sizeof(int) * n);
                    if (!xs || !ys) {
                        if (xs) free(xs); if (ys) free(ys);
                        if (rank == 0) fprintf(stderr, "Failed to allocate coords for n=%d\n", n);
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }
                }
            } else {
                allocated_coords = 1;
                xs = malloc(sizeof(int) * n);
                ys = malloc(sizeof(int) * n);
                if (!xs || !ys) {
                    if (xs) free(xs); if (ys) free(ys);
                    if (rank == 0) fprintf(stderr, "Failed to allocate coords for n=%d\n", n);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }

            /* populate random coordinates identically on all ranks when allocated */
            if (allocated_coords) {
                /* ensure same sequence across ranks */
                srand(seed);
                for (int i = 0; i < n; ++i) {
                    xs[i] = rand() % 101;
                    ys[i] = rand() % 101;
                }
                if (rank == 0) printf("Generated random coords for n=%d run=%d seed=%u\n", n, run, seed);
            }

            travelingSalesman ts;
            travelingSalesman_init(&ts, xs, ys, n);

            double seq_elapsed = -1.0;
            coordinateSet seqBest;
            seqBest.path = NULL;
            seqBest.distance = DBL_MAX;

            /* Run sequential brute-force only on rank 0 to get baseline timing when requested */
            if (rank == 0 && (size == 1 || force_seq)) {
                double s0 = MPI_Wtime();
                seqBest = travelingSalesman_SEQ_bruteForce(&ts);
                double s1 = MPI_Wtime();
                seq_elapsed = s1 - s0;
                printf("[n=%d run=%d] Sequential: time = %.6f s, best distance = %.6f\n", n, run, seq_elapsed, seqBest.distance);
                if (seqBest.path) { free(seqBest.path); seqBest.path = NULL; }
                /* free costMatrix and re-init to avoid stale data for parallel run */
                travelingSalesman_free(&ts);
                travelingSalesman_init(&ts, xs, ys, n);
            } else if (rank == 0 && size > 1 && !force_seq) {
                printf("[n=%d run=%d] Skipping sequential baseline (run with --seq to force)\n", n, run);
            }

            MPI_Barrier(MPI_COMM_WORLD);

            if (size > 1) {
                /* Run parallel brute-force on all ranks and time wall-clock (max across ranks) */
                double p0 = MPI_Wtime();
                coordinateSet parBest = travelingSalesman_PAR_bruteForce(&ts, rank, size);
                double p1 = MPI_Wtime();
                double local_par_time = p1 - p0;
                double par_elapsed = 0.0;
                MPI_Reduce(&local_par_time, &par_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

                if (rank == 0) {
                    double best_distance = (parBest.path) ? parBest.distance : DBL_MAX;
                    if (parBest.path) { printf("[n=%d run=%d] Parallel: time = %.6f s, best distance = %.6f\n", n, run, par_elapsed, parBest.distance); free(parBest.path); parBest.path = NULL; }
                    else { printf("[n=%d run=%d] Parallel: no path returned to root\n", n, run); }

                    /* write CSV line for this run (root only) */
                    if (csvf) {
                        /* build row piece by piece so missing values are empty fields */
                        fprintf(csvf, "%d,%d,%d,%u,", n, size, run, seed);
                        if (seq_elapsed > 0.0) fprintf(csvf, "%.6f,", seq_elapsed);
                        else fprintf(csvf, ",");
                        fprintf(csvf, "%.6f,", par_elapsed);
                        fprintf(csvf, "%.6f\n", best_distance);
                        fflush(csvf);
                    }
                }
            } else {
                /* size == 1: only sequential run was performed; do not run parallel.
                   Write CSV with seq_time filled and leave par_time empty. */
                if (rank == 0) {
                    double best_distance = seqBest.distance;
                    if (csvf) {
                        fprintf(csvf, "%d,%d,%d,%u,", n, size, run, seed);
                        if (seq_elapsed > 0.0) fprintf(csvf, "%.6f,", seq_elapsed);
                        else fprintf(csvf, ",");
                        /* leave par_time empty */
                        fprintf(csvf, ",%.6f\n", best_distance);
                        fflush(csvf);
                    }
                }
            }

            /* cleanup per-run */
            travelingSalesman_free(&ts);
            if (allocated_coords) { free(xs); free(ys); }

            MPI_Barrier(MPI_COMM_WORLD);
        } /* end runs */
    } /* end n loop */

    if (rank == 0 && csvf) fclose(csvf);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
