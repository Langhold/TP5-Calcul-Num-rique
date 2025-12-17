/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
	
	for(int i = 0; i<*kv; ++i){
		for(int j = 0; j<*lab; ++j){
			AB[i*(*la+*kv) + j] = 0;
		}
	}
	
	for (int i = 0; i<*la; ++i) {
		int indx =i* *lab+*kv;
		AB[indx] = -1;
		AB[indx+1] =  2;
		AB[indx+2] = -1;
	}
	AB[*kv] = 0;
	AB[*la**lab-1] = 0;
	
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
	for(int i = 0; i<*kv; ++i){
		for(int j = 0; j<*lab; ++j){
			AB[i*(*la+*kv) + j] = 0;
		}
	}
	
	for (int i = 0; i<*la; ++i) {
		int indx =i* *lab+*kv;
		AB[indx+1] = 1;
	}
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
    for (int i = 0; i < *la; ++i) {
        RHS[i] = 0.0;
    }
	RHS[0] += *BC0;
	RHS[*la-1] += *BC1;
}

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
	for (int i = 0; i<*la; ++i) {
		EX_SOL[i] = *BC0 + X[i] * (*BC1-*BC0);
	}
}

void set_grid_points_1D(double* x, int* la){
	const double h = 1/((double)*la+1);
	for (int i = 0; i<*la; ++i) {
		x[i] = (i+1)*h;
	}
}

double relative_forward_error(double* x, double* y, int* la){
    double diff_norm = 0.0;
    double x_norm = 0.0;

    for (int i = 0; i < *la; ++i) {
        const double diff = y[i] - x[i];
        diff_norm += diff * diff;
        x_norm += x[i] * x[i];
    }

    double ferr = sqrt(diff_norm / x_norm);
    return ferr;
}

int indexABCol(int i, int j, int *lab){
	const int row = i-j;
	int indx = j* *lab +row;
	return indx;
}

static inline int indexABColtridiag(int i, int j, int *lab, int *ku){
	const int row = *ku+i-j;
	int indx = j* *lab +row;       /* raw = ku+1+i-j, col = j, colmajor_indx = col* *lab + row*/
	return indx;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
	
	for(int k = 1; k<*n; ++k){
		int indx_kk = indexABColtridiag(k,k,lab,ku);
		if (AB[indx_kk] == 0.0f) {
			*info = -1;
			return *info;
		}
		
		for (int i = k+1; i<k+*kl; i++) {
			int indx_ik = indexABColtridiag(i,k,lab,ku);
			AB[indx_ik] *= 1/AB[indx_kk];
			for (int j =k+1; j<*n; j++) {
				int indx_ij = indexABColtridiag(i,j,lab,ku);
				int indx_jk = indexABColtridiag(j,k,lab,ku);
				AB[indx_ij] -= AB[indx_ik]  * AB[indx_jk];
			}
		}
	}
  return *info;
}

