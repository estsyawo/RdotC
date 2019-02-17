// multiply two matrices
/*
This file contains functions for matrix and vector operations and other useful
functions
*/

# include "matply.h"

/* xa - vectorised matrix A, 
  xb - vectorised matrix B,
  xab - product AB (vectorised)
  nra - number of rows in A
  nca - number of columns in B
  ncb - number of columns in B
*/

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
        v += pow(absval(x[i]), (double) *p);
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
