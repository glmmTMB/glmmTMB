namespace glmmtmb{
  template<class Type> Type dbetabinom(Type y, Type a, Type b, Type n, int give_log=0);
  template<class Type> Type dgenpois(Type y, Type theta, Type lambda, int give_log=0);
  double rtruncpois(int k, double mu);
  template<class Type> Type rgenpois(Type theta, Type lambda);
  template<class Type> Type rtruncated_genpois(Type theta, Type lambda);
  template<class Type> Type rtruncated_compois2(Type mean, Type nu);
  template<class Type> Type rtweedie(Type mu, Type phi, Type p);
  template<class Type> Type logit_invcloglog(Type x);
  template<class Type> Type logit_pnorm(Type x);


extern "C" {
    /* See 'R-API: entry points to C-code' (Writing R-extensions) */
	double Rf_logspace_sub (double logx, double logy);
	void   Rf_pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);
	SEXP compois_calc_var(SEXP mean, SEXP nu);
}

} // namespace glmmtmb

// Quantile functions needed to simulate from truncated distributions
// (soon to be obsolete ...)
double Rf_qnbinom(double p, double size, double prob, int lower_tail, int log_p);
double Rf_qpois(double p, double lambda, int lower_tail, int log_p);
	
