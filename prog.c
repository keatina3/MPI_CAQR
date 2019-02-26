#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include "utils.h"
#include "dense_orthogon.h"

int main(int argc, char *argv[]){
    int myid, nprocs;
    int rows, cols;
    DenseMatrix A, Q, R, W;
    srand48(1000);
    //srand48(time(NULL));
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); 
    
    rows = 1200, cols = 20;
    W = gen_mat(rows/4, cols, 1);
    
    if(myid==0){
        R = gen_mat(cols, cols, 0);
        Q = gen_mat(rows, cols, 0);
    }

    TSQR(&W, &Q, &R, nprocs, myid, MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
