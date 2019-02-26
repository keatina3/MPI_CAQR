#ifndef _UTILS_H_
#define _UTILS_H_

typedef struct
{
	int I;
	int J;
	double *vals;
	double **col_ptr;
} DenseMatrix, *DenseMatrix_arr;

DenseMatrix gen_mat(int I, int J, int rand);
void assign_rand(DenseMatrix *A);
void init_mat(DenseMatrix *A);
void init_upper_trian(DenseMatrix *A);
void free_mat(DenseMatrix *A);
void print_mat(DenseMatrix *A);

#endif
