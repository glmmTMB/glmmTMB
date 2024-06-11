#include <stdlib.h>
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" {

  SEXP compois_calc_var(SEXP mean, SEXP nu);
  SEXP omp_check();

  const static R_CallMethodDef R_CallDef[] = {
    TMB_CALLDEFS,
    CALLDEF(compois_calc_var, 2),
    CALLDEF(omp_check, 0),
    {NULL, NULL, 0}
  };

  void R_init_glmmTMB(DllInfo *dll)
  {
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
#ifdef TMB_CCALLABLES
    TMB_CCALLABLES("glmmTMB");
#endif
  }

}
