#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <mpi.h>

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

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--seq") == 0) force_seq = 1;
    }

    MPI_Init(&argc, &argv);                 /* Initialize the MPI environment */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   /* Get the rank of the process */
    MPI_Comm_size(MPI_COMM_WORLD, &size);   /* Get the number of processes */

    if (rank == 0) printf("Comparing sequential and parallel TSP on %d MPI ranks\n", size);

    /* city coordinates small set for algorithm testing */
    int small_n = 5;
    int small_xs[5] = {0,1,2,3,5};
    int small_ys[5] = {0,1,2,3,8};

    travelingSalesman ts;
    travelingSalesman_init(&ts, small_xs, small_ys, small_n);

    double seq_elapsed = -1.0;
    coordinateSet seqBest;
    seqBest.path = NULL;
    seqBest.distance = DBL_MAX;

    /* Run sequential brute-force only on rank 0 to get baseline timing.
       Skip this when running a multi-rank job unless --seq is provided to avoid OOM. */
    if (rank == 0 && (size == 1 || force_seq)) {
        double s0 = MPI_Wtime();
        seqBest = travelingSalesman_SEQ_bruteForce(&ts);
        double s1 = MPI_Wtime();
        seq_elapsed = s1 - s0;
        printf("Sequential: time = %.6f s, best distance = %.6f\n", seq_elapsed, seqBest.distance);
        if (seqBest.path) {
            printf("Sequential path: ");
            for (int i = 0; i < ts.numCities; ++i) {
                printf("%d", seqBest.path[i]);
                if (i + 1 < ts.numCities) printf(" -> ");
            }
            printf(" -> %d\n", seqBest.path[0]);
            free(seqBest.path);
            seqBest.path = NULL;
        }
        /* free costMatrix and any internal allocations so parallel run can rebuild/broadcast */
        travelingSalesman_free(&ts);
        travelingSalesman_init(&ts, small_xs, small_ys, small_n);
    } else if (rank == 0 && size > 1 && !force_seq) {
        printf("Skipping sequential baseline in multi-rank run. To force it, run with --seq (may OOM).\n");
    }

    /* synchronize all ranks before parallel run */
    MPI_Barrier(MPI_COMM_WORLD);

    /* Run parallel brute-force on all ranks and time wall-clock (max across ranks) */
    double p0 = MPI_Wtime();
    coordinateSet parBest = travelingSalesman_PAR_bruteForce(&ts, rank, size);
    double p1 = MPI_Wtime();
    double local_par_time = p1 - p0;
    double par_elapsed = 0.0;
    MPI_Reduce(&local_par_time, &par_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        if (parBest.path) {
            printf("Parallel: time = %.6f s, best distance = %.6f\n", par_elapsed, parBest.distance);
            printf("Parallel path: ");
            for (int i = 0; i < ts.numCities; ++i) {
                printf("%d", parBest.path[i]);
                if (i + 1 < ts.numCities) printf(" -> ");
            }
            printf(" -> %d\n", parBest.path[0]);
            free(parBest.path);
            parBest.path = NULL;
        } else {
            printf("Parallel: no path returned to root\n");
        }

        if (seq_elapsed > 0.0 && par_elapsed > 0.0) {
            double S = seq_elapsed / par_elapsed;
            double E = S / (double)size;
            printf("Speedup S = T_seq / T_par = %.4f\n", S);
            printf("Efficiency E = S / p = %.4f\n", E);
        } else {
            printf("Cannot compute speedup/efficiency (missing timings)\n");
        }
    }

    /* cleanup and exit */
    travelingSalesman_free(&ts);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
