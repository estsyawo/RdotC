#include <R.h>
#include <math.h>
#include "matply.h"
#define absval(x) ((x) >=0.0 ? (x):(-(x)))
// declare function prototypes.

// void solve(double *A, double *b, double *X, int *n, char *type);
void gauss_seidel( double *A, double *b, double *phi, int *n);
void conjgrad(double *A, double *b, double *X, int *n, double *r, double *p, double *vec);
void SOR( double *A, double *b, double *phi, int *n);
