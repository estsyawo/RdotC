#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void linreg_cg(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void linreg_cord(void *, void *, void *, void *, void *, void *, void *);
extern void linreg_gs(void *, void *, void *, void *, void *, void *, void *);
extern void linreg_qr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void linreg_sor(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"linreg_cg",   (DL_FUNC) &linreg_cg,   10},
  {"linreg_cord", (DL_FUNC) &linreg_cord,  7},
  {"linreg_gs",   (DL_FUNC) &linreg_gs,    7},
  {"linreg_qr",   (DL_FUNC) &linreg_qr,   13},
  {"linreg_sor",  (DL_FUNC) &linreg_sor,   7},
  {NULL, NULL, 0}
};

void R_init_RdotC(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
