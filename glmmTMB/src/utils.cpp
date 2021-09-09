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
    return ScalarLogical(1);
#else
    return ScalarLogical(0);
#endif
}

/* openmp controller */
extern "C"
SEXP omp_num_threads(SEXP x) {
#ifdef _OPENMP
  if( !isNull(x) ){
    int n = INTEGER(x)[0];
    omp_set_num_threads( n );
  }
  return ScalarInteger( omp_get_max_threads() );
#else
  warning("OpenMP not supported.");
  return ScalarInteger( 0 );
#endif
}
