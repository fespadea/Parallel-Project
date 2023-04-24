#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "clockcycle.h"
#include "proj.h"

#define clock_frequency 512000000

int main(int argc, char **argv){
    FILE* metaDataFile = fopen("formattedMetaData.dat", "r");
    int n;
    fscanf(metaDataFile, "%d", &n);
    int m;
    fscanf(metaDataFile, "%d", &m);
    fclose(metaDataFile);
    printf("%d %d\n", n, m);
    int totalLength = n * m;

    
    int rank, nprocs;
    MPI_File fh;
    MPI_Status status;
    int bufsize, nints;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    bufsize = totalLength*sizeof(int)/nprocs;
    int *buf = (int *)malloc(bufsize);
    nints = bufsize/sizeof(int);

    unsigned long long start_cycles= clock_now();
    
    printf("%i\n", buf[nints-1]);
    MPI_File_open(MPI_COMM_WORLD, "formattedData.dat", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_at(fh, rank*bufsize, buf, nints, MPI_INT, &status);
    MPI_File_close(&fh);

    printf("%i\n", buf[nints-1]);
    
    FILE* dataFile = fopen("formattedData.dat", "r");
    double ** dataMatrix = (double**)malloc(n * sizeof(double*));
    int localn = nints/m;
    for(int i = 0; i < localn; i++){
        dataMatrix[i] = (double*)malloc(m * sizeof(double));
        for(int j = 0; j < m; j++){
            dataMatrix[i][j] = buf[i*m + j];
        }
    }

    double epsilon = 0.05;
    double delta = 0.1;
    int sMult = 3;
    double alpha = 0.5;

    double ** ATilde = matrixSparsification(dataMatrix, localn, m, epsilon, delta, sMult, alpha, n, rank, nprocs);

    for(int i = 0; i < localn; i++){
        for(int j = 0; j < m; j++){
            buf[i*m + j] = ATilde[i][j];
        }
    }
    MPI_File_open(MPI_COMM_WORLD, "output.dat", MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_at(fh, rank*bufsize, buf, nints, MPI_INT, &status);
    MPI_File_close(&fh);

    unsigned long long end_cycles= clock_now();

    printf("%lf\n", error(dataMatrix, ATilde, n, m));
    double time_in_secs_CUDA = ((double)(start_cycles - end_cycles)) / clock_frequency;
    printf("CUDA Reduce Sum Seconds Taken: %lf\n", time_in_secs_CUDA);

    free(ATilde);
    free(dataMatrix);
    free(buf);


    MPI_Finalize();
}