/*
 Emmanuel S. Tsyawo
 estsyawo@temple.edu,  estsyawo@gmail.com
 February 08, 2019
 A wrapper function for solvers of system of linear equations
 */
/* Compile using
 gcc -c SOR_test.c matply.c
 gcc -o execSOR_test SOR_test.o matply.o
 ./execSOR_test
 */
//**************************************************************************************//
// This file pools together all solvers for systems of linear equations
//**************************************************************************************//

#include "solve.h"

/*
  Input:
	A - nxn matrix
	b - nx1 vector
	X - nx1 vector of unknowns to be solved for
	type - a character for type of solver. solvers currently available are
	"conjgrad", "gauss_seidel", and "SOR"
*/
//**************************************************************************************//
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

    // free allocated memory
    free(r);
    free(p);
    free(vec);
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
/*
 Compile using: SOR.c matply.c
 */


//**************************************************************************************//




//**************************************************************************************//
