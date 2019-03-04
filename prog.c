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
	double time_taken;
	clock_t start, end;
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

	A = gen_mat(rows, cols, 1);
	Q = gen_mat(rows, cols, 0);
	R = gen_mat(cols, cols, 0);

	start=clock();	
	CAQR(&A, &Q, &R, block, 0, A.I, R.I, nprocs, myid, MPI_COMM_WORLD);
	end=clock();
	time_taken = ((double)(end-start))/CLOCKS_PER_SEC;
	printf("%d, %d, %d, %lf\n",rows,cols,block,time_taken);
	
	free_mat(&A); free_mat(&Q); free_mat(&R);
	MPI_Finalize();

    return 0;
}
