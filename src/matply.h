/*header file for functions in matply.c*/
#include <R.h>
#include <math.h>
#define absval(x) ((x) >=0.0 ? (x):(-(x)))

// declare function prototypes.
void matply(double *xa, double *xb,double *xab, int *nra, int *nca,int *ncb);
void matply_sym(double *x, double *xx, int *nrx, int *ncx);
void matply_xty(double *x, double *y, double *xty, int *nrx, int *ncx, int *ncy);
void matply_xyt(double *x, double *y, double *xyt, int *nrx, int *ncx, int *nry);
void trans(double *a, double *atrans,int *nra,int *nca);
double dotprod(double *a, double *b, int *n);
void matply_ax(double *x, double *a, double *ax, int *n);
void vecadd(double *xa, double *xb, double *xm, int *n);
void vecsub(double *xa, double *xb, double *xm, int *n);
double max(double *x, int *n);
double min(double *x, int *n);
double norm_lp(double *x, int *n, int *p);
double norm_max(double *x, int *n);
