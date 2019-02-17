#include <R.h>
#include <math.h>
#define absval(x) ((x) >=0.0 ? (x):(-(x)))

// declare function prototypes.
void linreg_cg(double *Y, double *X, double *coefs, int *nrX, int *ncX, double *XX, double *XY);
void linreg_gs(double *Y, double *X, double *coefs, int *nrX, int *ncX, double *XX, double *XY);
void linreg_sor(double *Y, double *X, double *coefs, int *nrX, int *ncX, double *XX, double *XY);
