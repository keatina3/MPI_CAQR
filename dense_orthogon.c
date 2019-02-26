#include <stdio.h>

#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "utils.h"
#include "dense_orthogon.h"

void TSQR(DenseMatrix *W, DenseMatrix* Q, DenseMatrix *R, int nprocs, int myid, MPI_Comm comm){
    int qI = W->I, qJ = W->J;
    int i,j, count=0;
    MPI_Status stat;
    DenseMatrix_arr *w;
    DenseMatrix Q_loc1, Q_loc2, R_hat, R_upper;
    
    w = (DenseMatrix_arr*)malloc((1+log2(nprocs))*sizeof(DenseMatrix*));
    for(i=0;i<(1+log2(nprocs));i++){
        w[i] = (DenseMatrix*)malloc(sizeof(DenseMatrix));
    }
    //Q_loc1 = gen_mat(qI, qJ, 0);
    //Q_loc2 = gen_mat(2*qJ, qJ, 0);
    R_upper = gen_mat(qJ, qJ, 0);
    R_hat = gen_mat(2*qJ, qJ, 0);
    
    printf("TEST1\n");
    for(i=1; i<=nprocs; i*=2){
        if(myid%i == 0){
            if(i==1)
                hhorth(W, &Q_loc1, &R_upper, w[count]);
            else if(i==nprocs){
                hhorth(&R_hat, &Q_loc2, &R_upper, w[count]);
                break;
            } else
                hhorth(&R_hat, &Q_loc2, &R_upper, w[count]);
            sendR(&R_hat, &R_upper, i, myid, comm);
        }
        count++;
    }
    printf("TEST2\n");
    
    count=0; 
    for(i=1; i<=nprocs; i*=2){
        if(myid%i == 0)
            free_mat(w[count]);
        free(w[count]);
        count++;
    }
    free(w);
}

void sendR(DenseMatrix *R_hat, DenseMatrix *R_upper, int level, int myid, MPI_Comm comm){
    int j;
    MPI_Status stat;
    for(j=0;j<R_upper->J;j++){
        if(myid%(level*2) != 0)
            MPI_Send(R_upper->col_ptr[j], (j+1), MPI_DOUBLE, (myid-level), 0, comm);
        else
            MPI_Recv(&R_hat->col_ptr[j][R_upper->I], (j+1), 
                MPI_DOUBLE, (myid+level), MPI_ANY_TAG, comm, &stat);
    }
}


void hhorth(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R, DenseMatrix *W){
	int i,j;
	int n = A->I, m = A->J;
	double *w_vals = (double*)calloc(n*m,sizeof(double));
	double **w = (double**)malloc(n*sizeof(double*));		// to store matrix of w's
	double *z = (double*)calloc(n,sizeof(double));			// stores z vector
	double *wv = (double*)calloc(n,sizeof(double));			// temp storage for Wj*V)
	double v;
	double *Xk = (double*)calloc(n,sizeof(double));		// temp placeholder of X
	//DenseMatrix W;
	// initialising pointers for W //
	for(j=0;j<m;j++){
    	w[j] = &w_vals[j*n];
	}
    *W = gen_mat(n, m, 0);
	
    for(j=0;j<m;j++){
		for(i=0;i<n;i++)
			Xk[i] = A->col_ptr[j][i];			// setting Xj = Aj
        if(j > 0){
            for(i=0;i<j;i++){			// iterating through (Pj-1)(Pj-2)...(P1)(Xj)
                //v = 2.0*inner_prod(Xk, w[i], n);	// getting v=2*X(^T)w
                v = 2.0*inner_prod(Xk, W->col_ptr[i], n);	// getting v=2*X(^T)w
                //vec_scal_prod(wv, w[i], v, n, 0); 	// getting wv(^T)
                vec_scal_prod(wv, W->col_ptr[i], v, n, 0); 	// getting wv(^T)
                vec_vec_add(Xk, Xk, wv, n, 1);		// PXk = Xk - wv(^T)
            }
		}
		get_z(z, Xk, j, n);					// getting z to be used in getting Wj
        //vec_scal_prod(w[j], z, norm_2(z, n), n, 1);		// w = z / ||z||
        vec_scal_prod(W->col_ptr[j], z, norm_2(z, n), n, 1);		// w = z / ||z||
        //v = 2.0*inner_prod(Xk, w[j], n);		// see above
        v = 2.0*inner_prod(Xk, W->col_ptr[j], n);		// see above
        //vec_scal_prod(wv, w[j], v, n, 0);
        vec_scal_prod(wv, W->col_ptr[j], v, n, 0);
        vec_vec_add(Xk, Xk, wv, n, 1);
		for(i=0;i<=j;i++)
			R->col_ptr[j][i] = Xk[i];			//	rk = Pk(Pk-1)...Xk
        /*
        Q->col_ptr[j][j] = 1.0;					//	Qej
        for(i=j;i>=0;i--){						// Qk = P1..(Pk)(ek)
            v = 2.0*inner_prod(Q->col_ptr[j], w[i], n);
            vec_scal_prod(wv, w[i], v, n, 0); 
            vec_vec_add(Q->col_ptr[j], Q->col_ptr[j], wv, n, 1);
        }
        */
    }
	// freeing allocated memory //
	free(Xk);
    free(w_vals);
    free(w);
    free(z);
    free(wv);
}

// applies piecewise function of z //
void get_z(double *z, double *x,  int k, int n){
	int i,j;
	double beta, tmp = 0.0;

	for(i=0;i<n;i++){
		if(i<k){
			z[i] = 0;
		} else if(i==k){
			for(j=k;j<n;j++){
				tmp += x[j]*x[j];
			}
            beta = x[k] > 0.0 ? sqrt(tmp) : -sqrt(tmp);
			z[i] = beta + x[i];
		} else {
			z[i] = x[i];
		}
	}
}

// scalar product //
// division if div=TRUE //
void vec_scal_prod(double *xhat, double *x, double y, int n, int div){
	int i;
    if(fabs(y) < ERR)
		return;
	
	for(i=0;i<n;i++){
		if(div)
			xhat[i] = x[i]/y;
		else
			xhat[i] = x[i]*y;
	}
}

// vector-vector add //
// subtract if sub=TRUE //
void vec_vec_add(double *xhat, double *x, double *y, int n, int sub){
	int i;
	for(i=0;i<n;i++){
		if(sub)
			xhat[i] = x[i] - y[i];
		else 
			xhat[i] = x[i] + y[i];
	}
}


// matrix-matrix multiply //
void mat_mul(DenseMatrix *AB, DenseMatrix *A, DenseMatrix *B){
	int i,j,k;
	if(A->J != B->I)
		return;
	for(j=0;j<B->J;j++)
		for(i=0;i<A->I;i++)
			for(k=0;k<A->J;k++)
				AB->col_ptr[j][i] += A->col_ptr[k][i]*B->col_ptr[j][k];

}
/*
// matrix-matrix addition //
void mat_mat_add(DenseMatrix *A_B, DenseMatrix *A, DenseMatrix* B, int sub){
    int j;
    if(A->J != B->J || A->I != B->I)
        return;
    for(j=0;j<A_B->J;j++)
        vec_vec_add(A_B->col_ptr[j], A->col_ptr[j], B->col_ptr[j], A_B->I, sub);
}

// matrix transpose //
void transpose(DenseMatrix *A, DenseMatrix *At){
	int i,j;
	
	for(j=0;j<At->J;j++)
		for(i=0;i<At->I;i++)
			At->col_ptr[j][i] = A->col_ptr[i][j];
}
*/

// euclidian norm //
double norm_2(double *x, int n){
	double norm=0;
	int i;
	for(i=0;i<n;i++)
		norm += x[i]*x[i];
	norm = sqrt(norm);

	return norm;
}

// inner-product of 2 vectors //
double inner_prod(double *x, double *y, int n){
	double prod;
	int i;
	for(i=0;i<n;i++)
		prod += x[i]*y[i];
	
	return prod;
}

/*
// calculates max( |Aij| )
double fwd_err(DenseMatrix *A){
	int i,j;
	double fwd_err = 0.0;
	for(j=0;j<A->J;j++)
		for(i=0;i<A->I;i++)
			fwd_err = A->col_ptr[j][i] > fwd_err ? A->col_ptr[j][i] : fwd_err;
	return fwd_err;
}
*/
