#include "RdotC.h"
#include <math.h>

// multiply two matrices

// take a product xab = xa*xb
void matply(double *xa, double *xb,double *xab, int *nra, int *nca,int *ncb){
  double sum ;
  for(int i=0; i< *nra;i++){
    for (int j=0; j<*ncb; j++) {
      sum = 0.0;
      for (int k=0; k<*nca; k++) {
        sum =  (double) (sum + (xa[ k*(*nra)+i])*(xb[ j*(*nca)+k])) ;
      }
      xab[ j*(*nra)+i] = sum ;
    }
  }
}

// take product x'x using only x as input; output ncx x ncx
// Author: Clara-Christina Gerstner & Emmanuel S. Tsyawo
void matply_sym(double *x, double *xx, int *nrx, int *ncx)
{
  double sum ;
  int i,j,k;
  for( i=0; i< *ncx;i++){
    for ( j=0; j< *ncx; j++) {
      sum = 0.0;
      for ( k=0; k<*nrx; k++) {
        // take inner product of i and j columns of x
        sum += x[i*(*nrx)+k]*x[j*(*nrx)+k];
      }
      xx[ j*(*ncx)+i] = sum ;
      if (j!=i) {
        xx[i*(*ncx)+j]=xx[j*(*ncx)+i];
      }
    }
  }
}

/*
Take the product X'Y (dimension [ncx,ncy]) using input X,Y. Note: nrx = nry
Transpose is integrated in the function

// Author: Clara-Christina Gerstner & Emmanuel S. Tsyawo
*/
void matply_xty(double *x, double *y, double *xty, int *nrx, int *ncx, int *ncy)
{
  int i, j, k;
  double sum;
  for (i=0; i<*ncx; i++) {
    for (j=0; j<*ncy; j++) {
      sum=0.0;
      for (k=0; k<*nrx; k++) {
        sum += x[i*(*nrx)+k]*y[j*(*nrx)+k];
      }
      xty[j*(*ncx)+i] = sum;
    }
  }
}



/*
Take the product XY' (dimension [nrx,nry]) using input X,Y. Note: ncx = ncy
Transpose is integrated in the function

// Author: Clara-Christina Gerstner & Emmanuel S. Tsyawo
*/
void matply_xyt(double *x, double *y, double *xyt, int *nrx, int *ncx, int *nry)
{
  int i, j, k;
  double sum;
  for (i=0; i<*nrx; i++) {
    for (j=0; j<*nry; j++) {
      sum=0.0;
      for (k=0; k<*ncx; k++) {
        sum += x[k*(*nrx)+i]*y[k*(*nry)+j];
      }
      xyt[j*(*nrx)+i] = sum;
    }
  }
}

// take product xa*xb excluding column ica in xa and row ica in xb
void matply_sk1(double *xa, double *xb,double *xab, int *nra, int *nca,int *ncb, int *ica)
{
  double sum ;
  for(int i=0; i< *nra;i++){
    for (int j=0; j<*ncb; j++) {
      sum = 0.0;
      for (int k=0; k<*nca; k++) {
        if (k!=*ica) { //skip column ica in xa and row ica in xb
          sum =  (double) (sum + (xa[ k*(*nra)+i])*(xb[ j*(*nca)+k])) ;
        }
      }
      xab[ j*(*nra)+i] = sum ;
    }
  }
}

/*
Take the transpose of a matrix a, store in atrans,
*/
void trans(double *a, double *atrans, int *nra,int *nca)
{
  int i,j;
  for(i=0;i<*nra;i++){
    for(j=0;j<*nca;j++){
      atrans[i*(*nca) +j] = a[j*(*nra)+i];
    }
  }
}

// take the dot product two vectors a, b with length n each.
double dotprod(double *a, double *b, int *n)
{
  int i;
  double ans;
  ans = 0.0;
  for (i=0; i<*n; i++) {
    ans += a[i]*b[i];
  }
  return ans;
}

// compute dot products of two columns taken from matrices a and b of nr
// number of rows each. column indices are ica and icb
double dotprod_col_ex(double *a, double *b, int *nr, int *ica, int *icb)
{
  int i;
  double ans;
  ans = 0.0;
  for (i=0; i<*nr; i++) {
    ans += a[(*ica-1)*(*nr)+i]*b[(*icb-1)*(*nr)+i];
  }
  return ans;
}


// take the product of a scalar and matrix/vector ax
void matply_ax(double *x, double *a, double *ax, int *n)
{
  int i;
  for (i=0; i<*n; i++) {
    ax[i] = (*a)*x[i];
  }
}

// add matrices/vectors xm = xa + xb
void vecadd(double *xa, double *xb, double *xm, int *n)
{
  int i;
  for (i=0; i<*n; i++) {
    xm[i] = xa[i] + xb[i];
  }
}

// add matrices/vectors xm = xa - xb
void vecsub(double *xa, double *xb, double *xm, int *n)
{
  int i;
  for (i=0; i<*n; i++) {
    xm[i] = xa[i] - xb[i];
  }
}

// find the maximum in a vector x of length n
double max(double *x, int *n)
{
  double xmax=x[0]; // initialise with first element
  for(int i=1; i<*n; i++){
    if(x[i]>xmax){
      xmax=x[i];
    }
  }
  return xmax;
}

// find the minimum in a vector x of length n
double min(double *x, int *n)
{
  double xmin=x[0]; // initialise with first element
  for(int i=1; i<*n; i++){
    if(x[i]<xmin){
      xmin=x[i];
    }
  }
  return xmin;
}

// compute p-norm of a vector x
double norm_lp(double *x, int *n, int *p)
{
  double v=0.0;
  int i;
  for (i=0; i<*n; i++) {
    v += pow( absval(x[i]), (double) *p);
  }
  return pow(v,(double) 1/(*p));
}

// compute the infinity norm of a vector x
double norm_max(double *x, int *n)
{
  double v;
  double xmax=absval(x[0]); // initialise with first element
  for(int i=1; i<*n; i++){
    v = absval(x[i]);
    if(v>xmax){
      xmax=v;
    }
  }
  return xmax;
}

//**************************************************************************************//

// conjugate gradient method for solving a system of linear equations AX=b

// function for the conjugate gradient algorithm
void conjgrad(double *A, double *b, double *X, int *n, double *r, double *p, double *vec)
{
  double alf, beta, tol, dev, dpr1, dpr0;
  int ncX, i,k,maxiter;

  // set parameters
  ncX = 1; //number of columns in X
  tol = 1e-7;
  maxiter = 1000; // maximum number of iterations allowed.

  k = 0;
  for(i=0;i<*n;i++){ // fill in initial r
    X[i] = 0.0; // initialise X to zero => A*X = 0 => initial r=b
    r[i] = b[i];
    p[i] = r[i]; //initialise p with r
  }//end for i

  dpr0=dotprod(r, r, n); // take dot product
  // initialise main while loop
  for(;;){

    matply(A, p, vec, n, n, &ncX); // take matrix product vec=A*p
    dpr1=dotprod(p, vec, n);// take dot product, store in dpr1 temporarily

    alf = dpr0/dpr1;

    //update X
    for(i=0;i<*n;i++){ // update X, r
      X[i] += alf*p[i];
      r[i] += -alf*vec[i] ; // recall vec=Ap currently
      vec[i] = absval(r[i]); // pass |r| to vec
    }

    dev = max(vec, n); //compute infinity norm of r

    // checking for convergence
    if(dev<tol){
      break;
    }
    if(k>=maxiter){
      Rprintf("Maximum number of iterations reached. Conjugate gradient algorithm failed to converge.\n");
      break;
    }

    dpr1 = dotprod(r,r,n); // take dot product dpr1
    beta = dpr1/dpr0;

    for(i=0;i<*n;i++){// update p
      p[i] = r[i] + beta*p[i];
    }

    k += 1; //update number of iterations k
    dpr0=dpr1; // update drp0

  }// end for(;;)

}

/*Compile using the following files: conjgrad.c matply.c*/


//**************************************************************************************//

void gauss_seidel( double *A, double *b, double *x, int *n){
  double sig, vl, tol, dev, mxdev;
  int i, j, iter, maxiter;
  tol=1e-10;
  iter = 0;
  maxiter = 1000;

  for(;;) { //begin do while loop
    iter +=1;
    for(i=0; i<*n; i++){
      sig = 0.0; dev=0.0; mxdev=0.0;
      for(j=0; j<*n; j++ ){
        if(j!=i){
          sig = (double) sig + A[j*(*n) + i]*x[j];
        } // end if
      }// end for j

      if(absval(A[i*(1+(*n))])<tol){
        Rprintf("Algorithm stopped: non-dominant diagonal term");
        break;
      }
      vl = x[i];
      x[i] = (double) ((b[i]-sig) - x[i])/A[i*(1+(*n))];
      dev = absval((x[i]-vl));
      if(dev>mxdev){
        mxdev = dev;
      }
    }// end for i

    if(mxdev<=tol){
      break;
    }
    if (iter>=maxiter) {
      Rprintf("Maximum number of iterations reached. Algorithm failed to converge.");
      break;
    }

  }
}

//**************************************************************************************//

void SOR( double *A, double *b, double *phi, int *n){
  double sig, vl, tol, dev, mxdev, w;
  int i, j, iter, maxiter;
  tol=1e-7;
  w = 1.0; // this value can be adjusted on the interval (0,2)
  maxiter = 1000;
  iter = 0;
  for(;;) { //begin do while loop
    iter +=1;
    for(i=0; i<*n; i++){
      sig = 0.0;
      dev = 0.0; mxdev = 0.0; // reset in order to check convergence
      for(j=0; j<*n; j++ ){
        if(j!=i){
          sig = (double) sig + A[j*(*n) + i]*phi[j];
        } // end if
      }// end for j

      if(absval(A[i*(1+(*n))])<tol){
        break;
      }

      vl = (double) w*(((b[i]-sig)/A[i*(1+(*n))]) - phi[i]);
      phi[i] = phi[i] + vl;
      dev = absval(vl);
      if(dev>mxdev){
        mxdev = dev; // update infinity norm of deviations
      }
    }// end for i
    if(mxdev<=tol){ // check for convergence
      break;
    }
    if (iter>=maxiter) {
      Rprintf("Warning: Maximum number of iterations reached.\n");
      break;
    }

  }
}

//****************************************************************************//
void linreg_qr(double *x, int *n, int *p, double *y, int *ny,
             double *tol, double *b, double *rsd,
             double *qty, int *k,
             int *jpvt, double *qraux, double *work)
{
  F77_CALL(dqrls)(x, n, p, y, ny,tol, b, rsd, qty, k, jpvt, qraux, work);
}
//****************************************************************************//




//****************************************************************************//

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
//****************************************************************************//




//****************************************************************************//
// linear regression by coordinate descent
void linreg_cord(double *Y, double *X, double *coefs, int *nrX, int *ncX, double *Xdot, double *ervec)
{
  double az = 1e-20, tol = 1e-07,RSS0,RSS1,dev,dr;
  int ncY=1, k=0,j,id,maxk = 2000;

  // compute column-wise sum-of-squares
  for (j=2; j<=*ncX; j++) {
    Xdot[j-1]=dotprod_col_ex(X, X, nrX, &j, &j);
  }//end for
  Xdot[0] = (double) *nrX; //intercept of 1's sums squares to nrX
  // compute initial sum of squared of errors
  matply(X, coefs, ervec, nrX, ncX, &ncY); //compute ervec = X*coefs
  vecsub(Y, ervec, ervec, nrX); // compute ervec <-- Y-ervec
  RSS0=dotprod(ervec, ervec, nrX); // compute initial residual sum of squares (RSS0)

  // commence iteration
  for(;;)
  {
    k+=1; //increment counter by 1
    for (j=0; j<*ncX; j++){
      // compute X*coefs excluding column j in X and element j in coefs
      matply_sk1(X, coefs,ervec, nrX, ncX, &ncY, &j); // ervec = X[,-j]*coef[-j]
      vecsub(Y, ervec, ervec, nrX); // compute ervec <-- Y-ervec
      id = j+1; // dotprod_col_ex() needs id, that counts from 1,..,ncX
      dr=dotprod_col_ex(ervec, X, nrX, &ncY, &id); // compute dot(ervec,X[,j])
      coefs[j] = dr/Xdot[j]; // update coefs
    }// end for j

    // compute sum of squared residuals
    matply(X, coefs, ervec, nrX, ncX, &ncY); //compute ervec = X*coefs
    vecsub(Y, ervec, ervec, nrX); // compute ervec <-- Y-ervec
    RSS1=dotprod(ervec, ervec, nrX); // compute residual sum of squares (RSS1)

    // compute deviation
    dev = (RSS0-RSS1)/(az + RSS0);

    // check for convergence
    if (dev<=tol) {
      break;
    }
    if (k>=maxk) {
      Rprintf("Warning: maximum number of iterations reached. \n");
      break;
    }
    RSS0 = RSS1; //update RSS0
  }// end for(;;)

}

