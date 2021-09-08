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

/* Swap x$obj and y$obj (assumes obj is 1st component!) */
extern "C"
SEXP obj_swap(SEXP x, SEXP y) {
  SEXP x_obj = VECTOR_ELT(x, 0);
  SEXP y_obj = VECTOR_ELT(y, 0);
  SET_VECTOR_ELT(x, 0, y_obj);
  SET_VECTOR_ELT(y, 0, x_obj);
  return R_NilValue;
}
