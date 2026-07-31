// additional distributions etc. for glmmTMB
#include <vector> // for vectors in Bell()

namespace glmmtmb{
	
  /* Not used anymore: */
  template<class Type>
  Type dbetabinom(Type y, Type a, Type b, Type n, int give_log=0)
  {
    /*
      Wikipedia:
      f(k|n,\alpha,\beta) =
      \frac{\Gamma(n+1)}{\Gamma(k+1)\Gamma(n-k+1)}
      \frac{\Gamma(k+\alpha)\Gamma(n-k+\beta)}{\Gamma(n+\alpha+\beta)}
      \frac{\Gamma(\alpha+\beta)}{\Gamma(\alpha)\Gamma(\beta)}
    */
    Type logres =
      lgamma(n + 1) - lgamma(y + 1)     - lgamma(n - y + 1) +
      lgamma(y + a) + lgamma(n - y + b) - lgamma(n + a + b) +
      lgamma(a + b) - lgamma(a)         - lgamma(b) ;
    if(!give_log) return exp(logres);
    else return logres;
  }
  /* R:
    > identical(lgamma(exp(-150)), 150)
    [1] TRUE
    FIXME: Move 'logspace_gamma' to TMB.
  */
  namespace adaptive {
    template<class T>
    T logspace_gamma(const T &x) {
      /* Tradeoff: The smaller x the better approximation *but* the higher
         risk of psigamma() overflow */
      if (x < -150)
        return -x;
      else
        return lgamma(exp(x));
    }
  }
  TMB_BIND_ATOMIC(logspace_gamma, 1, adaptive::logspace_gamma(x[0]))
  template<class Type>
  Type logspace_gamma(Type x) {
    CppAD::vector<Type> args(2); // Last index reserved for derivative order
    args[0] = x;
    args[1] = 0;
    return logspace_gamma(args)[0];
  }
  template<class Type>
  Type dbetabinom_robust(Type y, Type loga, Type logb, Type n, int give_log=0)
  {
    Type a = exp(loga), b = exp(logb);
    Type logy = log(y), lognmy = log(n - y); // May be -Inf
    Type logres =
      lgamma(n + 1) - lgamma(y + 1) - lgamma(n - y + 1) +
      logspace_gamma(logspace_add(logy, loga)) +
      logspace_gamma(logspace_add(lognmy, logb)) -
      lgamma(n + a + b) +
      lgamma(a + b) - logspace_gamma(loga) - logspace_gamma(logb);
    if(!give_log) return exp(logres);
    else return logres;
  }
	
  template<class Type>
  Type dgenpois(Type y, Type theta, Type lambda, int give_log=0)
  {
    /*
      f(y|\theta,\lambda) =
      \frac{\theta(theta+\lambda y)^{y-1}e^{-\theta-\lambda y}}{y \!}
    */
    Type logres =
      log(theta) + (y - 1) * log(theta + lambda * y) -
      theta - lambda * y - lgamma(y + Type(1));
    if(!give_log) return exp(logres);
    else return logres;
  }

  template<class Type>
  Type dskewnorm(Type y, Type mu, Type sigma, Type alpha, int give_log=0)
  {
    Type delta = alpha/sqrt(1 + pow(alpha, 2)); 
    Type omega = sigma/sqrt(1 - 2/M_PI * pow(delta, 2)); 
    Type xi = mu - omega * delta * sqrt(2/M_PI); 
    Type logres = 
      log(2.0) - log(omega) + log(dnorm((y - xi)/omega, Type(0), Type(1), 0)) + log(pnorm(alpha * (y - xi)/omega));
    if(!give_log) return exp(logres);
    else return logres;
  }
  
   // from C. Geyer aster package, src/raster.c l. 175
   // Simulate from truncated poisson
   // see https://cran.r-project.org/web/packages/aster/vignettes/trunc.pdf for technical/mathematical details
   // k is the truncation point (e.g. k=0 -> 0-truncated)
   // MODIFICATIONS: change die() to throw std::range_error() 
   double rtruncated_poisson(int k, double mu)
   {
    int m;
    double mdoub;

    if (mu <= 0.0)
	    throw std::range_error("non-positive mu in k-truncated-poisson simulator\n");
    if (k < 0)
	    throw std::range_error("negative k in k-truncated-poisson simulator\n");

    mdoub = k + 1 - mu;
    if (mdoub < 0.0)
        mdoub = 0.0;
    m = mdoub;
    if (m < mdoub)
        m = m + 1;
    /* since mu > 0.0 we have 0.0 <= mdoub < k + 1 hence 0 <= m <= k + 1 */

    for (;;) {
        double x = rpois(mu) + m;
        if (m > 0) {
            double a = 1.0;
            int j;
            double u = unif_rand();
            for (j = 0; j < m; ++j)
                a *= (k + 1 - j) / (x - j);
            if (u < a && x > k)
                return x;
        } else {
            if (x > k)
                return x;
        }
    }
   } // rtruncpois

   // alpha = size (dispersion param), k = truncation point, mu = mean
   double rtruncated_nbinom(double alpha, int k, double mu)
   {
    int m;
    double mdoub;
    double p = alpha / (mu + alpha);
    double q = mu / (mu + alpha);

    if (alpha <= 0.0)
        throw std::range_error("non-positive size in k-truncated-neg-bin simulator\n");
    if (mu <= 0.0)
        throw std::range_error("non-positive mu in k-truncated-neg-bin simulator\n");
    if (k < 0)
        throw std::range_error("negative k in k-truncated-neg-bin simulator\n");

    mdoub = (k + 1.0) * p - alpha * q;
    if (mdoub < 0.0)
        mdoub = 0.0;
    m = mdoub;
    if (m < mdoub)
        m = m + 1;
    /* since p < 1.0 and q > 0.0 we have 0.0 <= mdoub < k + 1
       hence 0 <= m <= k + 1 */

    for (;;) {
        double x = rnbinom(alpha + m, p) + m;
        if (m > 0) {
            double a = 1.0;
            int j;
            double u = unif_rand();
            for (j = 0; j < m; ++j)
                a *= (k + 1 - j) / (x - j);
            if (u < a && x > k)
                return x;
        } else {
            if (x > k)
                return x;
        }
    }
   } // rtruncated_nbinom

  /* Simulate from generalized poisson distribution */
  template<class Type>
  Type rgenpois(Type theta, Type lambda) {
    // Copied from R function HMMpa::rgenpois
    Type ans = Type(0);
    Type random_number = runif(Type(0), Type(1));
    Type kum = dgenpois(Type(0), theta, lambda);
    while (random_number > kum) {
      ans = ans + Type(1);
      kum += dgenpois(ans, theta, lambda);
    }
    return ans;
  }
	
  /* Simulate from zero-truncated generalized poisson distribution */
  template<class Type>
  Type rtruncated_genpois(Type theta, Type lambda) {
    int nloop = 10000;
    int counter = 0;
    Type ans = rgenpois(theta, lambda);
    while(ans < Type(1) && counter < nloop) {
      ans = rgenpois(theta, lambda);
      counter++;
    }
    if(ans < 1.) Rf_warning("Zeros in simulation of zero-truncated data. Possibly due to low estimated mean.");
    return ans;
  }

  extern "C" {
    /* See 'R-API: entry points to C-code' (Writing R-extensions) */
    double Rf_logspace_sub (double logx, double logy);
    void   Rf_pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);
  }
  /* y(x) = logit_invcloglog(x) := log( exp(exp(x)) - 1 ) = logspace_sub( exp(x), 0 )

     y'(x) = exp(x) + exp(x-y) = exp( logspace_add(x, x-y) )

   */
  TMB_ATOMIC_VECTOR_FUNCTION(
                             // ATOMIC_NAME
                             logit_invcloglog
                             ,
                             // OUTPUT_DIM
                             1,
                             // ATOMIC_DOUBLE
                             ty[0] = Rf_logspace_sub(exp(tx[0]), 0.);
                             ,
                             // ATOMIC_REVERSE
                             px[0] = exp( logspace_add(tx[0], tx[0]-ty[0]) ) * py[0];
                             )
  template<class Type>
  Type logit_invcloglog(Type x) {
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return logit_invcloglog(tx)[0];
  }

  /* y(x) = logit_pnorm(x) := logit( pnorm(x) ) =
     pnorm(x, lower.tail=TRUE,  log.p=TRUE) -
     pnorm(x, lower.tail=FALSE, log.p=TRUE)

     y'(x) = dnorm(x) * ( (1+exp(y)) + (1+exp(-y)) )

  */
  double logit_pnorm(double x) {
    double log_p_lower, log_p_upper;
    Rf_pnorm_both(x, &log_p_lower, &log_p_upper, 2 /* both tails */, 1 /* log_p */);
    return log_p_lower - log_p_upper;
  }
  TMB_ATOMIC_VECTOR_FUNCTION(
                             // ATOMIC_NAME
                             logit_pnorm
                             ,
                             // OUTPUT_DIM
                             1,
                             // ATOMIC_DOUBLE
                             ty[0] = logit_pnorm(tx[0])
                             ,
                             // ATOMIC_REVERSE
                             Type zero = 0;
                             Type tmp1 = logspace_add(zero, ty[0]);
                             Type tmp2 = logspace_add(zero, -ty[0]);
                             Type tmp3 = logspace_add(tmp1, tmp2);
                             Type tmp4 = dnorm(tx[0], Type(0), Type(1), true) + tmp3;
                             px[0] = exp( tmp4 ) * py[0];
                             )
  template<class Type>
  Type logit_pnorm(Type x) {
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return logit_pnorm(tx)[0];
  }

  /* Calculate variance in compois family using

     V(X) = (logZ)''(loglambda)

  */
  double compois_calc_var(double mean, double nu){
    using atomic::compois_utils::calc_loglambda;
    using atomic::compois_utils::calc_logZ;
    double loglambda = calc_loglambda(log(mean), nu);
    typedef atomic::tiny_ad::variable<2, 1, double> ADdouble;
    ADdouble loglambda_ (loglambda, 0);
    ADdouble ans = calc_logZ<ADdouble>(loglambda_, nu);
    return ans.getDeriv()[0];
  }

  /* Calculate variance in combinom family using

     V(X) = (logZ)''(Psi)

  */
  double combinom_calc_var(double mean, double nu, int n){
    using atomic::combinom_utils::calc_logitp;
    using atomic::combinom_utils::calc_logZ;
    // NB calc_logitp takes the mean itself, not log(mean)
    // (unlike compois_utils::calc_loglambda above)
    double logitp = calc_logitp(mean, nu, n);
    typedef atomic::tiny_ad::variable<2, 1, double> ADdouble;
    ADdouble logitp_ (logitp, 0);
    ADdouble ans = calc_logZ<ADdouble>(logitp_, nu, n);
    return ans.getDeriv()[0];
  }

  /* Simulate from zero-truncated Conway-Maxwell-Poisson distribution */
  template<class Type>
  Type rtruncated_compois2(Type mean, Type nu) {
    int nloop = 10000;
    int counter = 0;
    Type ans = rcompois2(mean, nu);
    while(ans < 1. && counter < nloop) {
      ans = rcompois2(mean, nu);
      counter++;
    }
    if(ans < 1.) Rf_warning("Zeros in simulation of zero-truncated data. Possibly due to low estimated mean.");
    return ans;
  }

  /* Simulate from Conway-Maxwell-Binomial distribution.
     CMB has finite support {0,...,n}. We build the unnormalized log-PMF
     once using the adjacent-ratio recurrence
       log w_k = log w_{k-1} + theta + nu * log((n-k+1)/k)
     then sample via inverse CDF. This avoids repeated dcombinom2 calls
     (each of which would recompute Z internally) and reduces simulation
     from O(n^2) to O(n) per draw. */
  template<class Type>
  Type rcombinom2(Type mean, Type nu, Type size) {
    int n = (int) asDouble(size);
    if (n <= 0) return Type(0);

    double theta = asDouble(combinom_calc_logitp(mean, nu, size));
    double nu_d  = asDouble(nu);

    // Unnormalized log-weights via adjacent-ratio recurrence
    std::vector<double> logw(n + 1);
    logw[0] = 0.0;
    double max_logw = 0.0;
    for (int k = 1; k <= n; k++) {
      logw[k] = logw[k - 1] + theta +
                nu_d * std::log((double)(n - k + 1) / (double)k);
      if (logw[k] > max_logw) max_logw = logw[k];
    }

    // Stabilized unnormalized total
    double total = 0.0;
    for (int k = 0; k <= n; k++) total += std::exp(logw[k] - max_logw);

    // Inverse-CDF draw, no PMF normalization needed
    double u = asDouble(runif(Type(0), Type(1))) * total;
    double cum = 0.0;
    for (int k = 0; k <= n; k++) {
      cum += std::exp(logw[k] - max_logw);
      if (cum >= u) return Type(k);
    }
    return Type(n);
  }

  /* Conway-Maxwell-Binomial PMF and CDF (plain double, mean-parameterized)
     backing the R-level dcombinom()/pcombinom(). Builds the unnormalized
     log-PMF over the finite support {0,...,n} with the same adjacent-ratio
     recurrence as rcombinom2, then normalizes. Returns NaN outside the
     valid parameter range (0 < mean < n). */
  double combinom_pmf_cdf(double x, double mean, double nu, int n,
                          bool cumulative, bool give_log) {
    if (n <= 0 || !R_FINITE(mean) || !R_FINITE(nu) ||
        mean <= 0. || mean >= (double) n) return R_NaN;
    double xf = std::floor(x);
    if (cumulative) {
      if (!(xf >= 0.)) return give_log ? R_NegInf : 0.;  // includes NaN x
      if (xf >= (double) n) return give_log ? 0. : 1.;
    } else {
      if (!(x >= 0.) || x > (double) n || x != xf)
        return give_log ? R_NegInf : 0.;
    }
    using atomic::combinom_utils::calc_logitp;
    double theta = calc_logitp(mean, nu, n);
    std::vector<double> logw(n + 1);
    logw[0] = 0.0;
    double max_logw = 0.0;
    for (int k = 1; k <= n; k++) {
      logw[k] = logw[k - 1] + theta +
                nu * std::log((double)(n - k + 1) / (double)k);
      if (logw[k] > max_logw) max_logw = logw[k];
    }
    double total = 0.0;
    for (int k = 0; k <= n; k++) total += std::exp(logw[k] - max_logw);
    double logtotal = max_logw + std::log(total);
    if (!cumulative) {
      double ans = logw[(int) xf] - logtotal;
      return give_log ? ans : std::exp(ans);
    }
    double cum = 0.0;
    for (int k = 0; k <= (int) xf; k++) cum += std::exp(logw[k] - logtotal);
    double p = std::min(std::max(cum, 0.), 1.);
    return give_log ? std::log(p) : p;
  }

  /* Simulate from tweedie distribution */
  template<class Type>
  Type rtweedie(Type mu, Type phi, Type p) {
    // Copied from R function tweedie::rtweedie
    Type lambda = pow(mu, 2. - p) / (phi * (2. - p));
    Type alpha  = (2. - p) / (1. - p);
    Type gam = phi * (p - 1.) * pow(mu, p - 1.);
    int N = (int) asDouble(rpois(lambda));
    Type ans = rgamma(N, -alpha /* shape */, gam /* scale */).sum();
    return ans;
  }

  /* Simulate from skew-normal distribution */
  template<class Type>
  Type rskewnorm(Type mu, Type sigma, Type alpha) {
    // Copied from R function sn::rsn
    Type delta = alpha/sqrt(1 + pow(alpha, 2)); 
    Type omega = sigma/sqrt(1 - 2/M_PI * pow(delta, 2)); 
    Type xi = mu - omega * delta * sqrt(2/M_PI); 
    
    Type chi = CppAD::abs(rnorm(Type(0), Type(1)));
    Type nrv = rnorm(Type(0), Type(1));
    Type z = delta * chi + sqrt(1 - pow(delta, 2)) * nrv;
    Type ans = xi + omega * z;
    
    return ans;
  }

  // Simulate from Bell distribution
  // translated from bellreg::rbell
  template<class Type>
  Type rbell(Type theta) {

    Type ans = 0;
    double dtheta = asDouble(theta);
    
    double lambda = expm1(dtheta);
    int N = (int) asDouble(rpois(lambda));
    for (int i=0; i<N; i++) {
      ans +=  rtruncated_poisson(0, dtheta);
    }
    
    return ans;
  }


  // FIXME: check!
  template<class Type>
  Type dcauchy(Type x, Type loc, Type scale, int give_log=0)
 {
   Type resid = (x - loc) / scale;
   Type logans = Type(-log(M_PI)) - log(scale) - log1p(resid*resid);
   if(give_log) return logans; else return exp(logans);
 }
  VECTORIZE4_ttti(dcauchy)

// first few Bell numbers, from https://oeis.org/A000110: 1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975, 678570, 4213597, 27644437, 190899322, 1382958545, 10480142147, 82864869804, 682076806159, 5832742205057, 51724158235372, 474869816156751, 4506715738447323, 44152005855084346, 445958869294805289, 4638590332229999353, 49631246523618756274
// implement as lookup table?

  double Bell(int n) {
  
    if ((n == 0) || (n == 1)) {
      return(1);
    }
    vector<double> B(n);
    vector<double> Bneu(n);
    B[0] = 1;
    for (int i=0; i<(n-1); i++) {
      Bneu[0] = B[i];
      for (int j=1; j<=(i+1); j++) {
	Bneu[j] = B[j - 1] + Bneu[j - 1];
      }
      for (int k=0; k<n; k++) B[k] = Bneu[k];
    }
    return Bneu[n-1];
  }

  template<class Type>
  Type dbell(Type x, Type theta, int give_log = 0)
  {

    // TMB doesn't have expm1 ... could use
    // exp(logspace_sub(theta, Type(0))) but feels like overkill,
    // maybe not what we want anyway?
    Type logres = x * log(theta) - exp(theta) + 1 +
      // clunky cast (ADvar -> double -> int)
      log(Bell((int) asDouble(x))) - lgamma(x + 1);
    if (!give_log) return exp(logres);
    return logres;
  }

  // taken from TMB atomic functions example,
  //    https://kaskr.github.io/adcomp/AtomicFunctions.html
  // if we wanted/needed to speed this up we could make use of the implementation
  // of Fukushima 2013 doi:10.1016/j.cam.2012.11.021 from 
  //     https://github.com/DarkoVeberic/LambertW

  // Double version of Lambert W function
double LambertW(double x) {
  double logx = log(x);
  double y = (logx > 0 ? logx : 0);
  int niter = 100, i=0;
  for (; i < niter; i++) {
    if ( fabs( logx - log(y) - y) < 1e-9) break;
    y -= (y - exp(logx - y)) / (1 + y);
  }
  if (i == niter) Rf_warning("W: failed convergence");
  return y;
}
  
TMB_ATOMIC_VECTOR_FUNCTION(
    // ATOMIC_NAME
    LambertW
    ,
    // OUTPUT_DIM
    1,
    // ATOMIC_DOUBLE
    ty[0] = LambertW(tx[0]); // Call the 'double' version
    ,
    // ATOMIC_REVERSE
    Type W  = ty[0];                    // Function value from forward pass
    Type DW = 1. / (exp(W) * (1. + W)); // Derivative
    px[0] = DW * py[0];                 // Reverse mode chain rule
)

// Scalar version
template<class Type>
Type LambertW(Type x){
  CppAD::vector<Type> tx(1);
  tx[0] = x;
  return LambertW(tx)[0];
}
// Vectorized version
VECTORIZE1_t(LambertW)


} // namespace glmmtmb

/* Interface to compois variance */
extern "C" {
  SEXP compois_calc_var(SEXP mean, SEXP nu) {
    if (LENGTH(mean) != LENGTH(nu))
      error("'mean' and 'nu' must be vectors of same length.");
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, LENGTH(mean)));
    for(int i=0; i<LENGTH(mean); i++)
      REAL(ans)[i] = glmmtmb::compois_calc_var(REAL(mean)[i], REAL(nu)[i]);
    UNPROTECT(1);
    return ans;
  }

  SEXP combinom_calc_var(SEXP mean, SEXP nu, SEXP size) {
    if (LENGTH(mean) != LENGTH(nu) || LENGTH(mean) != LENGTH(size))
      error("'mean', 'nu', and 'size' must be vectors of same length.");
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, LENGTH(mean)));
    for(int i=0; i<LENGTH(mean); i++)
      REAL(ans)[i] = glmmtmb::combinom_calc_var(
        REAL(mean)[i], REAL(nu)[i], INTEGER(size)[i]);
    UNPROTECT(1);
    return ans;
  }

  /* all arguments must already be recycled to a common length in R */
  SEXP dcombinom_R(SEXP x_, SEXP mean_, SEXP nu_, SEXP size_, SEXP give_log_) {
    int n = LENGTH(x_);
    if (LENGTH(mean_) != n || LENGTH(nu_) != n || LENGTH(size_) != n)
      error("'x', 'mean', 'nu', and 'size' must be vectors of same length.");
    int give_log = asLogical(give_log_);
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, n));
    for (int i = 0; i < n; i++)
      REAL(ans)[i] = glmmtmb::combinom_pmf_cdf(
        REAL(x_)[i], REAL(mean_)[i], REAL(nu_)[i], INTEGER(size_)[i],
        false, give_log);
    UNPROTECT(1);
    return ans;
  }

  SEXP pcombinom_R(SEXP q_, SEXP mean_, SEXP nu_, SEXP size_) {
    int n = LENGTH(q_);
    if (LENGTH(mean_) != n || LENGTH(nu_) != n || LENGTH(size_) != n)
      error("'q', 'mean', 'nu', and 'size' must be vectors of same length.");
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, n));
    for (int i = 0; i < n; i++)
      REAL(ans)[i] = glmmtmb::combinom_pmf_cdf(
        REAL(q_)[i], REAL(mean_)[i], REAL(nu_)[i], INTEGER(size_)[i],
        true, false);
    UNPROTECT(1);
    return ans;
  }
}

/* Interface to generalized Poisson sampler, Bell PMF, and Lambert W */
extern "C" {
  /* lambda1_ and lambda2_ must already be recycled to length n in R */
  SEXP rgenpois_R(SEXP n_, SEXP lambda1_, SEXP lambda2_) {
    int n = asInteger(n_);
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, n));
    GetRNGstate();
    for (int i = 0; i < n; i++)
      REAL(ans)[i] = glmmtmb::rgenpois(REAL(lambda1_)[i], REAL(lambda2_)[i]);
    PutRNGstate();
    UNPROTECT(1);
    return ans;
  }

  SEXP dbell_R(SEXP x_, SEXP theta_, SEXP give_log_) {
    int n = LENGTH(x_);
    if (LENGTH(theta_) != n)
      error("'x' and 'theta' must have the same length");
    int give_log = asLogical(give_log_);
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, n));
    for (int i = 0; i < n; i++)
      REAL(ans)[i] = glmmtmb::dbell(REAL(x_)[i], REAL(theta_)[i], give_log);
    UNPROTECT(1);
    return ans;
  }

  SEXP lambertW_R(SEXP x_) {
    int n = LENGTH(x_);
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, n));
    for (int i = 0; i < n; i++)
      REAL(ans)[i] = glmmtmb::LambertW(REAL(x_)[i]);
    UNPROTECT(1);
    return ans;
  }
}


