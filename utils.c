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
	
	init_mat(&A,A.I);
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
void init_mat(DenseMatrix *A, int offset){
	int j;
	for(j=0;j<A->J;j++)
		A->col_ptr[j] = &A->vals[(offset)*j];
}

int decomp1d(int n, int p, int myid, int *s, int *e){
    int d,r;
    d = n/p;
    r = n%p;
    if(myid < r){
    	*s = myid*(d+1);
    	*e = *s + d;
    } else {
    	*s = r*(d+1) + (myid-r)*d;
    	*e = *s + (d-1);
    }
 	return 0;
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
