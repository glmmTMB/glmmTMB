#include <stdlib.h>
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" {

  SEXP compois_calc_var(SEXP mean, SEXP nu);
  SEXP combinom_calc_var(SEXP mean, SEXP nu, SEXP size);
  SEXP omp_check();
  SEXP rgenpois_R(SEXP n_, SEXP lambda1_, SEXP lambda2_);
  SEXP dbell_R(SEXP x_, SEXP theta_, SEXP give_log_);
  SEXP lambertW_R(SEXP x_);
  SEXP dcombinom_R(SEXP x_, SEXP mean_, SEXP nu_, SEXP size_, SEXP give_log_);
  SEXP pcombinom_R(SEXP q_, SEXP mean_, SEXP nu_, SEXP size_);

  const static R_CallMethodDef R_CallDef[] = {
    TMB_CALLDEFS,
    CALLDEF(compois_calc_var, 2),
    CALLDEF(combinom_calc_var, 3),
    CALLDEF(omp_check, 0),
    CALLDEF(rgenpois_R, 3),
    CALLDEF(dbell_R, 3),
    CALLDEF(lambertW_R, 1),
    CALLDEF(dcombinom_R, 5),
    CALLDEF(pcombinom_R, 4),
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
