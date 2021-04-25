# include <R.h>
# include <Rinternals.h>
# ifdef _OPENMP
#include <omp.h>
# endif

/* check if openmp is enabled */
extern "C"
SEXP omp_check(void) {
#ifdef _OPENMP
    return ScalarLogical(1);
#else
    return ScalarLogical(0);
#endif
}
