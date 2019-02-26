#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <mpi.h>
#include <time.h>
#include "utils.h"
#include "dense_orthogon.h"

int main(int argc, char *argv[]){
    int myid, nprocs;
    int option;
    int rows, cols, block, t;
    DenseMatrix A, Q, R, W;
    srand48(1000);
    //srand48(time(NULL));
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); 
    
    rows = 1E5, cols = 5E4, block = 500, t=0;
	
    while((option=getopt(argc,argv,"n:m:b:t"))!=-1){
		switch(option){
			case 'n': rows = atoi(optarg);
				break;
			case 'm': cols = atoi(optarg);
				break;
			case 'b': block = atoi(optarg);
				break;
			case 't': t = 1;
				break;
			default:
				printf("Incorrect options entered!\n");
				return 1;
		}
	}	
	if(argc != optind){
		printf("Too many arguments provided, exiting!\n");
		return 1;
	}

    W = gen_mat(rows/4, block, 1);

    if(myid==0){
        R = gen_mat(cols, cols, 0);
        Q = gen_mat(rows, cols, 0);
    }

    TSQR(&W, &Q, &R, nprocs, myid, MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
