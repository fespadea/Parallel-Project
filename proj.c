#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <mpi.h>
#include "clockcycle.h"
#include "proj.h"

#define clock_frequency 512000000

int main(int argc, char **argv){
    FILE* metaDataFile = fopen("formattedMetaData.tsv", "r");
    int n;
    fscanf(metaDataFile, "%d", &n);
    int m;
    fscanf(metaDataFile, "%d", &m);
    fclose(metaDataFile);
    printf("%d %d\n", n, m);
    int totalLength = n * m;

    unsigned long long start_cycles= clock_now();
    
    FILE* dataFile = fopen("formattedData.tsv", "r");
    double ** dataMatrix = (double**)malloc(n * sizeof(double*));
    for(int i = 0; i < n; i++){
        dataMatrix[i] = (double*)malloc(m * sizeof(double));
        for(int j = 0; j < m; j++){
            fscanf(dataFile, "%lf", &dataMatrix[i][j]);
        }
    }

    double epsilon = 0.05;
    double delta = 0.1;
    int sMult = 3;
    double alpha = 0.5;

    double ** ATilde = matrixSparsification(dataMatrix, n, m, epsilon, delta, sMult, alpha);

    unsigned long long end_cycles= clock_now();

    printf("%lf\n", error(dataMatrix, ATilde, n, m));

    free(ATilde);
    free(dataMatrix);

    
    // // MPI init stuff
    // MPI_Init(&argc, &argv);
    // int myrank;
    // int numranks;
    // MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    // MPI_Comm_size(MPI_COMM_WORLD, &numranks);


    

    // MPI_Finalize();
}