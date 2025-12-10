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
		int indx =i* *lab+*kv +1;
		AB[indx] = 1;
		printf("ah !");
	}
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // TODO: Compute RHS vector
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // TODO: Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D
}  

void set_grid_points_1D(double* x, int* la){
  // TODO: Generate uniformly spaced grid points in [0,1]
}

double relative_forward_error(double* x, double* y, int* la){
  // TODO: Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)
  return 0.0;
}

int indexABCol(int i, int j, int *lab){
  // TODO: Return the correct index formula for column-major band storage
  return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // TODO: Implement specialized LU factorization for tridiagonal matrices
  return *info;
}
