#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include "utils.h"
#include "dense_orthogon.h"

int main(int argc, char *argv[]){
    int myid, nprocs;
    srand48(1000);
    //srand48(time(NULL));
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    

    MPI_Finalize();

    return 0;
}
