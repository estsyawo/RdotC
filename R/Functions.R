# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#===============================================================================#
# Notes :
# when Eror in FUN(X[[i]],...) : shows, run devtools::load_all(".")
#
# 1. use the current R Development Version (that will eventually become 3.4)
# 2. Run the tools::package_native_routine_registration_skeleton(".") and copy and
#    paste the full output in a packagename_init.c file to be put in src/
# 3. update NAMESPACE, verifying that useDynLib(packagename, .registration = TRUE)
# 4. If necessary, replace the exportPattern with export( list of object to be exported )
#===============================================================================#

#' Linear regression with Conjugate Gradient Solver
#'
#' This function implements linear regression by solving the normal equations. It
#' constructs matrices X'X and X'Y and solves for beta in X'X*beta=X'Y
#' using the iterative Conjugate Gradient solver of system of linear equations.
#'
#' @param Y vector of outcome variable
#' @param X matrix of covariates
#' @return coefs vector of coefficients
#'
#' @examples
#' linreg.cg(women$height,women$weight)
#' @useDynLib RdotC linreg_cg
#' @export

linreg.cg<- function(Y,X){
  if(is.numeric(Y)&&is.numeric(X)){
  X = as.matrix(cbind(1,X)) # add intercept term
  ncX = as.integer(ncol(X)); nrX = as.integer(nrow(X))
  ans=.C("linreg_cg",as.double(Y),as.double(X),coefs=double(ncX),nrX,ncX,
         XX=double(as.integer(ncX^2)),double(ncX),double(ncX),double(ncX),double(ncX),
         PACKAGE = "RdotC")
  # linreg_cg(double *Y, double *X, double *coefs, int *nrX, int *ncX, double *XX,
  # double *XY, double *r, double *p, double *vec)
  list(coefs=ans$coefs)
  }else{
    stop("Data input must be numeric.")
  }
}

#===============================================================================#

#' Linear regression with the SOR Solver
#'
#' This function implements linear regression by solving the normal equations. It
#' constructs matrices X'X and X'Y and solves for beta in X'X*beta=X'Y
#' using the iterative SOR solver of system of linear equations.
#'
#' @param Y vector of outcome variable
#' @param X matrix of covariates
#' @return coefs vector of coefficients
#'
#' @examples
#' linreg.sor(women$height,women$weight)
#' @useDynLib RdotC linreg_sor
#' @export

linreg.sor<- function(Y,X){
  if(is.numeric(Y)&&is.numeric(X)){
    X = as.matrix(cbind(1,X)) # add intercept term
    ncX = as.integer(ncol(X)); nrX = as.integer(nrow(X))
    ans=.C("linreg_sor",as.double(Y),as.double(X),coefs=double(ncX),nrX,ncX,
           XX=double(as.integer(ncX^2)),XY=double(ncX),
           PACKAGE = "RdotC")
    # linreg_cg(double *Y, double *X, double *coefs, int *nrX, int *ncX, double *XX,
    # double *XY, double *r, double *p, double *vec)
    list(coefs=ans$coefs)
  }else{
    stop("Data input must be numeric.")
  }
}

#===============================================================================#

#' Linear regression with the Gauss-Seidel Solver
#'
#' This function implements linear regression by solving the normal equations. It
#' constructs matrices X'X and X'Y and solves for beta in X'X*beta=X'Y
#' using the iterative Gauss-Seidel solver of system of linear equations.
#'
#' @param Y vector of outcome variable
#' @param X matrix of covariates
#' @return coefs vector of coefficients
#'
#' @examples
#' linreg.sor(women$height,women$weight)
#' @useDynLib RdotC linreg_gs
#' @export

linreg.gs<- function(Y,X){
  if(is.numeric(Y)&&is.numeric(X)){
    X = as.matrix(cbind(1,X)) # add intercept term
    ncX = as.integer(ncol(X)); nrX = as.integer(nrow(X))
    ans=.C("linreg_gs",as.double(Y),as.double(X),coefs=double(ncX),nrX,ncX,
           XX=double(as.integer(ncX^2)),XY=double(ncX),
           PACKAGE = "RdotC")
    # linreg_cg(double *Y, double *X, double *coefs, int *nrX, int *ncX, double *XX,
    # double *XY, double *r, double *p, double *vec)
    list(coefs=ans$coefs)
  }else{
    stop("Data input must be numeric.")
  }
}

#===============================================================================#

#' Linear regression with Lapack's QR solver \code{dqrls}
#'
#' This function implements linear regression using QR decomposition.
#' Lapack fortran subroutine \code{dqrls} is called in C for this purpose.
#'
#' @param Y vector of outcome variable
#' @param X matrix of covariates
#' @return coefs vector of coefficients
#'
#' @examples
#' linreg.qr(women$height,women$weight)
#' @useDynLib RdotC linreg_qr
#' @export

linreg.qr<- function(Y,X)
{
  X = as.matrix(cbind(1,X)) # add intercept term
  n = as.integer(nrow(X)); p = as.integer(ncol(X)); ny = 1; tol=1e-7
  X = as.double(X); Y = as.double(Y); ny = as.integer(ny); tol = as.double(tol)
  ans=.C("linreg_qr",X,n,p,Y,ny,tol,coefs=double(p),residuals=double(n),double(n),
         integer(ny),jpvt=integer(p),double(p),double(as.integer(2*p)))
  list(coefs=ans$coefs,residuals=ans$residuals,jpvt=ans$jpvt)
}

#===============================================================================#

#' Linear regression with Coordinate descent
#'
#' This function implements linear regression using the coordinate descent
#' algorithm. This minimises the sum of squares with respect to each parameter
#' at a time holding others fixed.
#'
#' @param Y vector of outcome variable
#' @param X matrix of covariates
#' @return coefs vector of coefficients
#'
#' @examples
#' linreg.coord(women$height,women$weight)
#' @useDynLib RdotC linreg_cord
#' @export

linreg.coord<- function(Y,X){
  if(is.numeric(Y)&&is.numeric(X)){
    X = as.matrix(cbind(1,X)) # add intercept term
    ncX = as.integer(ncol(X)); nrX = as.integer(nrow(X))
    ans=.C("linreg_cord",as.double(Y),as.double(X),coefs=double(ncX),nrX,ncX,
           Xdot=double(ncX),double(nrX),PACKAGE = "RdotC")
    list(coefs=ans$coefs)
  }else{
    stop("Data input must be numeric.")
  }
}

#===============================================================================#
