#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

// function to allocate and assign vals to new dense matrix //
// defaults to 0.0 if rand = FALSE //
DenseMatrix gen_mat(int I, int J, int rand){
	DenseMatrix A;
	
	A.I = I;
	A.J = J;
	
	A.vals = (double*)calloc(A.I*A.J,sizeof(double));
	A.col_ptr = (double**)malloc(A.J*sizeof(double*));
	
	init_mat(&A);
	if(rand)
		assign_rand(&A);

	return A;
}

// assigns random vals to matrix A //
void assign_rand(DenseMatrix *A){
	int i;
	for(i=0;i<A->I*A->J;i++)
		A->vals[i] = drand48();
}

// initialises pointers in 2d array in col major format //
void init_mat(DenseMatrix *A){
	int j;
	for(j=0;j<A->J;j++)
		A->col_ptr[j] = &A->vals[(A->I)*j];
}

// free allocated memory of dense matrix //
void free_mat(DenseMatrix *A){
	free(A->vals);
	free(A->col_ptr);
}

// printing matrix //
void print_mat(DenseMatrix *A){
	int i,j;
	for(i=0;i<A->I;i++){
		for(j=0;j<A->J;j++){
			printf("%0.6lf  ",A->col_ptr[j][i]);
		}
		printf("\n");
	}
}
