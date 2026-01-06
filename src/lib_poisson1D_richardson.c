/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

// TODO: Compute all eigenvalues for the 1D Poisson operator
void eig_poisson1D(double* eigval, int *la){
	double h = 1.0f/((double)*la+1.0f);
	for (int i = 0; i<*la; ++i) {
		eigval[i] = 2-2*cos(i*M_PI*h);
	}
}

double eigmax_poisson1D(int *la){
  return 4;
}

double eigmin_poisson1D(int *la){
  return 0;
}

double richardson_alpha_opt(int *la){
	const double temp = eigmax_poisson1D(la) + eigmin_poisson1D(la);
	double alpha = 2/temp;
	return alpha;
}

/**
 * Solve linear system Ax=b using Richardson iteration with fixed relaxation parameter alpha.
 * The iteration is: x^(k+1) = x^(k) + alpha*(b - A*x^(k))
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    const int n = *la;
    *nbite = 0;
	
    double *r = (double*)malloc((size_t)n * sizeof(double));
	if (!r){
		printf("Error: Allocation Failed\n");
		exit(1);
	}
	
    for (int i = 0; i < n; ++i) r[i] = RHS[i];
	cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1, AB, *lab, X, 1, 1, r, 1);

    double normb = cblas_dnrm2(*la, RHS, 1);
    int k = 0;
	
    while (k < *maxit) {
		
        double ferr = cblas_dnrm2(*la, r, 1) / normb;
        resvec[k] = ferr;
        if (ferr <= *tol) break;

        cblas_daxpy(*la, *alpha_rich, r, 1, X, 1);
        for (int i = 0; i < n; ++i) r[i] = RHS[i];
		cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku,-1.0, AB, *lab, X, 1, 1.0, r, 1);
        ++k;
		
    }
    *nbite = k;
    free(r);
}


/**
 * Extract MB for Jacobi method from tridiagonal matrix.
 * Such as the Jacobi iterative process is: x^(k+1) = x^(k) + D^(-1)*(b - A*x^(k))
 */
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
    for (int i = 0; i < *la; ++i) {
        MB[i] = AB[indexABColtridiag(i, i, lab, ku)];
    }
}

/**
 * Extract MB for Gauss-Seidel method from tridiagonal matrix.
 * Such as the Gauss-Seidel iterative process is: x^(k+1) = x^(k) + (D-E)^(-1)*(b - A*x^(k))
 */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
	
	int ii, jj, kk;
	for (jj=0;jj<(*la);jj++){
		kk = jj*(*lab);
		if (*kv>=0){
			for (ii=0;ii< *kv;ii++){
				MB[kk+ii]=0.0;
			}
		}
		MB[kk+ *kv]=0;
		MB[kk+ *kv+1]=AB[kk+ *kv+1];
		MB[kk+ *kv+2]=-AB[kk+ *kv+2];
	}
	if (*kv == 1) {MB[1]=0;}
	
	MB[(*lab)*(*la)-1]=0.0;
}

/**
 * Solve linear system Ax=b using preconditioned Richardson iteration.
 * The iteration is: x^(k+1) = x^(k) + M^(-1)*(b - A*x^(k))
 * where M is either D for Jacobi or (D-E) for Gauss-Seidel.
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
	const int n = *la;
	*nbite = 0;
	
	int * ipiv = (int *) calloc(*la, sizeof(int));
	
	double *r = (double*)malloc((size_t)n * sizeof(double));
	double *z = (double*)malloc((size_t)n * sizeof(double));
	double *temp = (double*)malloc((size_t)(*lab*n) * sizeof(double));
	if (!r){
		printf("Error: Allocation Failed\n");
		exit(1);
	}
	
	for (int i = 0; i < n; ++i) r[i] = RHS[i];
	cblas_dgbmv(CblasColMajor, CblasNoTrans, n, n, *kl, *ku, -1, AB, *lab, X, 1, 1, r, 1);

	double normb = cblas_dnrm2(n, RHS, 1);
	int k = 0;
	
	while (k < *maxit) {
		
		double ferr = cblas_dnrm2(n, r, 1) / normb;
		resvec[k] = ferr;
		if (ferr <= *tol) break;
		
		
		for (int i = 0; i < *lab*n; ++i) temp[i] = MB[i];
		
		LAPACKE_dgbsv(LAPACK_COL_MAJOR, *la, *kl, *ku, 1, MB, *lab, ipiv, r, *la); //M-1*r=z
		
		for (int i = 0; i < *lab; ++i) MB[i] = temp[i];
		for (int i = 0; i < n; ++i) z[i] = r[i];
		
		cblas_daxpy(*la, 1.0, z, 1, X, 1); //x = x + z
		for (int i = 0; i < n; ++i) r[i] = RHS[i];
		cblas_dgbmv(CblasColMajor, CblasNoTrans, n, n, *kl, *ku,-1.0, AB, *lab, X, 1, 1.0, r, 1);
		++k;
		
	}
	*nbite = k;
	free(r);
	free(z);
	free(temp);
	free(ipiv);
}

