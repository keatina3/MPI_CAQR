#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include "utils.h"
#include "dense_orthogon.h"

void CAQR(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R, int b, int count, int offsetA, int offsetR, int nprocs, int myid, MPI_Comm comm){
	MPI_Status stat;
	DenseMatrix Wi, R_upper, A_trail;
	DenseMatrix_arr *w;
	double *Rvals, *Avals;
	double **Rcols, **Acols;
	int i, j;
	int s, e, diff;
	int count2=0;
	
	if(count==0){
		Rvals = R->vals;
		Rcols = R->col_ptr;
		Avals = A->vals;
		Acols = A->col_ptr;
	}
	
	decomp1d(A->I, nprocs, myid, &s, &e);
	Wi.I = e-s+1;
	Wi.J = b;
	Wi.col_ptr = (double**)malloc(Wi.J*sizeof(double*));
	Wi.vals = &(A->vals[s]);
	init_mat(&Wi, offsetA);
	
	R_upper.I = b;
	R_upper.J = b;
	R_upper.col_ptr = (double**)malloc(R_upper.J*sizeof(double*));
	R_upper.vals = R->vals;
	init_mat(&R_upper, offsetR);		

	MPI_Barrier(comm);
	TSQR(&Wi, Q, &R_upper, &w, nprocs, myid, comm);
	
	A_trail.I = e-s+1;
	A_trail.J = (A->J)-b;
	A_trail.col_ptr = (double**)malloc(A_trail.J*sizeof(double*));
	A_trail.vals = &(A->col_ptr[b][s]);
	init_mat(&A_trail, offsetA);
	
	updateTrailing(&A_trail, &w, nprocs, myid);	
   
	for(i=1; i<=nprocs; i*=2){
        if(myid%i == 0)
            free_mat(w[count2]);
        free(w[count2]);
        count2++;
    }
	free(w);
	
	if(myid==0){
		for(i=0;i<A_trail.J;i++)
			memcpy(R->col_ptr[i+b], A_trail.col_ptr[i], b*sizeof(double));
	}
	
	for(i=1;i<nprocs;i++)
		if(myid==i)
			for(j=0;j<A_trail.J;j++)
				MPI_Send(A_trail.col_ptr[j], b, MPI_DOUBLE, (i-1), 0, comm);  
	for(i=0;i<(nprocs-1);i++)
		if(myid==i)
			for(j=0;j<A_trail.J;j++)
				MPI_Recv(&A_trail.col_ptr[j][A_trail.I], b, MPI_DOUBLE, (i+1), 0, comm, &stat);

	free(A_trail.col_ptr);
	free(Wi.col_ptr);
	free(R_upper.col_ptr);
	
	A->I = (A->I)-b;
	A->J = (A->J)-b;
	A->vals = &(A->col_ptr[b][b]);
	init_mat(A, offsetA);
	
	R->I = (R->I)-b;
	R->J = (R->J)-b;
	R->vals = &(R->col_ptr[b][b]);
	init_mat(R, offsetR);
	
	MPI_Barrier(comm);
	count++;
	if(myid==0)
		diff = (e-s+1) > b ? b : (e-s+1);
	MPI_Bcast(&diff, 1, MPI_INT, 0, comm);
	if(diff > nprocs && A->J > b){
		CAQR(A, Q, R, diff, count, offsetA, offsetR, nprocs, myid, comm);
	} else {
		if(myid==0)
			hhorth(A, Q, R, Q);
	}
	
	R->vals = Rvals;
	R->col_ptr = Rcols;
	A->vals = Avals;
	A->col_ptr = Acols;
}

void updateTrailing(DenseMatrix *A_tilde, DenseMatrix_arr **w, int nprocs, int myid){
	int i,j,k;
	double v, *wv;
	int n = A_tilde->I, m = A_tilde->J, count=0; 
	wv = (double*)calloc(n,sizeof(double));

	for(i=1;i<=nprocs;i*=2){
		if(myid%i == 0){
			for(j=0;j<m;j++){
    			for(k=0;k<(*w)[count]->J;k++){
					v = 2.0*inner_prod(A_tilde->col_ptr[j], (*w)[count]->col_ptr[k], n);
        			vec_scal_prod(wv, (*w)[count]->col_ptr[k], v, n, 0); 
            		vec_vec_add(A_tilde->col_ptr[j], A_tilde->col_ptr[j], wv, n, 1);
				}
			}
    	}
		count++;
	}
	free(wv);
}

void TSQR(DenseMatrix *W, DenseMatrix* Q, DenseMatrix *R, DenseMatrix_arr **w, int nprocs, int myid, MPI_Comm comm){
    int qJ = W->J;
    int i, count=0;
    DenseMatrix Q_loc1, Q_loc2, R_hat;
    
    (*w) = (DenseMatrix_arr*)malloc((1+log2(nprocs))*sizeof(DenseMatrix*));
    for(i=0;i<(1+log2(nprocs));i++){
        (*w)[i] = (DenseMatrix*)malloc(sizeof(DenseMatrix));
    }
    //Q_loc1 = gen_mat(qI, qJ, 0);
    //Q_loc2 = gen_mat(2*qJ, qJ, 0);
    //R_upper = gen_mat(qJ, qJ, 0);
	R_hat = gen_mat(2*qJ, qJ, 0);
    
    for(i=1; i<=nprocs; i*=2){
        if(myid%i == 0){
            if(i==1){
                hhorth(W, &Q_loc1, R, (*w)[count]);
            }else if(i==nprocs){
                hhorth(&R_hat, &Q_loc2, R, (*w)[count]);
                break;
            } else
                hhorth(&R_hat, &Q_loc2, R, (*w)[count]);
            sendR(&R_hat, R, i, myid, comm);
        }
        count++;
    }
	free_mat(&R_hat);
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
	double *z = (double*)calloc(n,sizeof(double));			// stores z vector
	double *wv = (double*)calloc(n,sizeof(double));			// temp storage for Wj*V)
	double v;
	double *Xk = (double*)calloc(n,sizeof(double));		// temp placeholder of X
    
	*W = gen_mat(n, m, 0);
	for(j=0;j<m;j++){
		for(i=0;i<n;i++)
			Xk[i] = A->col_ptr[j][i];			// setting Xj = Aj
		if(j > 0){
            for(i=0;i<j;i++){			// iterating through (Pj-1)(Pj-2)...(P1)(Xj)
                v = 2.0*inner_prod(Xk, W->col_ptr[i], n);	// getting v=2*X(^T)w
                vec_scal_prod(wv, W->col_ptr[i], v, n, 0); 	// getting wv(^T)
                vec_vec_add(Xk, Xk, wv, n, 1);		// PXk = Xk - wv(^T)
            }
		}
		get_z(z, Xk, j, n);					// getting z to be used in getting Wj
        vec_scal_prod(W->col_ptr[j], z, norm_2(z, n), n, 1);		// w = z / ||z||
        v = 2.0*inner_prod(Xk, W->col_ptr[j], n);		// see above
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
