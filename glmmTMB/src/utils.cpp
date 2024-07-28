/* as per BDR/CRAN, system headers must come before R headers */
# ifdef _OPENMP
#include <omp.h>
# endif
# include <R.h>
# include <Rinternals.h>

/* check if openmp is enabled */
extern "C"
SEXP omp_check(void) {
#ifdef _OPENMP
    return Rf_ScalarLogical(1);
#else
    return Rf_ScalarLogical(0);
#endif
}
