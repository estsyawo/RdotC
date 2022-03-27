#include <R.h>
#include <Rinternals.h>
#include <math.h>
#define absval(x) ((x) >=0.0 ? (x):(-(x)))
#define signum(x) ((x) >=0.0 ? (1):(-(1)))

// declare function prototypes.
void linreg_cg(double *Y, double *X, double *coefs, int *nrX, int *ncX, double *XX, double *XY, double *r, double *p, double *vec);
void linreg_gs(double *Y, double *X, double *coefs, int *nrX, int *ncX, double *XX, double *XY);
void linreg_sor(double *Y, double *X, double *coefs, int *nrX, int *ncX, double *XX, double *XY);
void F77_NAME(dqrls)(double *x, int *n, int *p, double *y, int *ny,
              double *tol, double *b, double *rsd,
              double *qty, int *k,
              int *jpvt, double *qraux, double *work);
// a c wrapper to F77_CALL(dqrls) - solve a least squares problem using QR decomposition
void linreg_qr(double *x, int *n, int *p, double *y, int *ny,
             double *tol, double *b, double *rsd,
             double *qty, int *k,
             int *jpvt, double *qraux, double *work);


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

void gauss_seidel( double *A, double *b, double *phi, int *n);
void conjgrad(double *A, double *b, double *X, int *n, double *r, double *p, double *vec);
void SOR( double *A, double *b, double *phi, int *n);
void Kern_Esc(double *z, double *Omg, int *n, int *ncz);
void Kern_DL(double *z, double *Omg, int *n, int *ncz);
