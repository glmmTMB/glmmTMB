# include <R.h>
# include <Rinternals.h>
# ifdef _OPENMP
#include <omp.h>
# endif

/* check if openmp is enabled */
SEXP ompcheck(void) {
#ifdef _OPENMP
    return ScalarLogical(1);
#else
    return ScalarLogical(0);
#endif
}

