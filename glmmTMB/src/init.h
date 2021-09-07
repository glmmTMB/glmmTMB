#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* FIXME: Won't be needed in upcoming TMB versions */
#ifndef TMB_CALLDEFS
#define TMB_CALLDEFS                                            \
  {"MakeADFunObject",     (DL_FUNC) &MakeADFunObject,     4},   \
  {"InfoADFunObject",     (DL_FUNC) &InfoADFunObject,     1},   \
  {"EvalADFunObject",     (DL_FUNC) &EvalADFunObject,     3},   \
  {"MakeDoubleFunObject", (DL_FUNC) &MakeDoubleFunObject, 3},   \
  {"EvalDoubleFunObject", (DL_FUNC) &EvalDoubleFunObject, 3},   \
  {"getParameterOrder",   (DL_FUNC) &getParameterOrder,   3},   \
  {"MakeADGradObject",    (DL_FUNC) &MakeADGradObject,    3},   \
  {"MakeADHessObject2",   (DL_FUNC) &MakeADHessObject2,   4},   \
  {"usingAtomics",        (DL_FUNC) &usingAtomics,        0},   \
  {"TMBconfig",           (DL_FUNC) &TMBconfig,           2}
#endif

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" {

  SEXP compois_calc_var(SEXP mean, SEXP nu);
  SEXP omp_check();
  SEXP omp_num_threads(SEXP);

  const static R_CallMethodDef R_CallDef[] = {
    TMB_CALLDEFS,
    CALLDEF(compois_calc_var, 2),
    CALLDEF(omp_check, 0),
    CALLDEF(omp_num_threads, 1),
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
