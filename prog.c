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
    DenseMatrix_arr *ws;
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

	/*	
    W = gen_mat(rows/4, block, 1);
	R = gen_mat(block, block, 0);
    Q = gen_mat(rows, cols, 0);

    TSQR(&W, &Q, &R, &ws, nprocs, myid, MPI_COMM_WORLD);

	free_mat(&W); free_mat(&Q); free_mat(&R);
    */
	
	A = gen_mat(rows, cols, 1);
	Q = gen_mat(rows, cols, 0);
	R = gen_mat(cols, cols, 0);
	
	CAQR(&A, &Q, &R, block, 0, A.I, R.I, nprocs, myid, MPI_COMM_WORLD);

	free_mat(&A); free_mat(&Q); free_mat(&R);
	MPI_Finalize();

    return 0;
}
