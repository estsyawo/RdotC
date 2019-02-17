/* 
Run linear regression for beta in Y = X*coefs + e. 
Solvers of system of normal equations available
 are 1. conjugate gradient, 2. gauss-seidel, and 3. SOR
*/

#include "matply.h"
#include "solve.h"

// linear regression with the conjugate gradient solver
void linreg_cg(double *Y, double *X, double *coefs, int *nrX, int *ncX, double *XX, double *XY, double *r, double *p, double *vec)
{
    int ncY = 1;
    
    matply_sym(X, XX, nrX, ncX); // take X'X
    matply_xty(X, Y, XY, nrX, ncX, &ncY); // take X'Y
    conjgrad(XX, XY, coefs, ncX, r, p, vec); // solve X'X*beta = X'Y
}

// linear regression with the gauss-seidel solver
void linreg_gs(double *Y, double *X, double *coefs, int *nrX, int *ncX, double *XX, double *XY)
{
    int ncY = 1;
    
    matply_sym(X, XX, nrX, ncX); // take X'X
    matply_xty(X, Y, XY, nrX, ncX, &ncY); // take X'Y
    gauss_seidel( XX, XY, coefs, ncX); // solve X'X*beta = X'Y
}

// linear regression with the sor solver
void linreg_sor(double *Y, double *X, double *coefs, int *nrX, int *ncX, double *XX, double *XY)
{
    int ncY = 1;
    
    matply_sym(X, XX, nrX, ncX); // take X'X
    matply_xty(X, Y, XY, nrX, ncX, &ncY); // take X'Y
    SOR(XX, XY, coefs, ncX); // solve X'X*beta = X'Y
}
