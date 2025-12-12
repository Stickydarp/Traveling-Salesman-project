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
uses DFS to find the shortest path
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

    MPI_Init(&argc, &argv);                 /* Initialize the MPI environment */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   /* Get the rank of the process */
    MPI_Comm_size(MPI_COMM_WORLD, &size);   /* Get the number of processes */

    printf("Hello from process %d of %d\n", rank, size);

    /* city coordinates small set for algorithm testing */
    int small_n = 5;
    int small_xs[5] = {0,1,2,3,5};
    int small_ys[5] = {0,1,2,3,8};

    travelingSalesman ts;
    travelingSalesman_init(&ts, small_xs, small_ys, small_n);

    double t0 = MPI_Wtime();
    coordinateSet trueBest = travelingSalesman_SEQ_bruteForce(&ts);
    double t1 = MPI_Wtime();

    printf("---- Rank %d/%d results ----\n", rank + 1, size);
    printf("Cities (%d):\n", ts.numCities);
    for (int i = 0; i < ts.numCities; ++i) {
        printf("  %2d: (%d, %d)\n", i, ts.xcoords[i], ts.ycoords[i]);
    }

    if (trueBest.path) {
        printf("Best distance: %.6f\n", trueBest.distance);
        printf("Path: ");
        for (int i = 0; i < ts.numCities; ++i) {
            printf("%d", trueBest.path[i]);
            if (i + 1 < ts.numCities) printf(" -> ");
        }
        if (ts.numCities > 0) printf(" -> %d", trueBest.path[0]); /* return to start */
        printf("\n");

        free(trueBest.path); /* caller frees returned path */
    } else {
        printf("No path found.\n");
    }

    printf("Compute time (this rank): %.6f s\n", t1 - t0);

    travelingSalesman_free(&ts);  /* free internal allocations */

    MPI_Barrier(MPI_COMM_WORLD);  /* synchronize before exit */
    MPI_Finalize();               /* Finalize the MPI environment */
    return 0;
}
