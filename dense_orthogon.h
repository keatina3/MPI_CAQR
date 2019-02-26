#ifndef _DENSE_ORTHOGON_H_
#define _DENSE_ORTHOGON_H_

#define ERR 1.0E-08

void TSQR(DenseMatrix *W, DenseMatrix* Q, DenseMatrix *R, int nprocs, int myid, MPI_Comm comm);
void sendR(DenseMatrix *R_hat, DenseMatrix *R_upper, int level, int myid, MPI_Comm comm);
void hhorth(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R);
void get_z(double *z, double *x, int k, int n);
void vec_scal_prod(double *x, double *xhat, double y, int n, int div);
void vec_vec_add(double *x, double *xhat, double *y, int n, int sub);
//void mat_mul(DenseMatrix *AB, DenseMatrix *A, DenseMatrix *B);
//void mat_mat_add(DenseMatrix *A_B, DenseMatrix *A, DenseMatrix *B, int sub);
//void transpose(DenseMatrix *A, DenseMatrix *At);
double norm_2(double *x, int n);
double inner_prod(double *x, double *y, int n);
//double fwd_err(DenseMatrix* A);

#endif
