#include <TMB.hpp>
#include "init.h"

namespace glmmtmb{
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
    if(ans < 1.) warning("Zeros in simulation of zero-truncated data. Possibly due to low estimated mean.");
    return ans;
  }

  template<class Type>
  bool isNA(Type x){
    return R_IsNA(asDouble(x));
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
    if(ans < 1.) warning("Zeros in simulation of zero-truncated data. Possibly due to low estimated mean.");
    return ans;
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
}

/* Interface to compois variance */
extern "C" {
  SEXP compois_calc_var(SEXP mean, SEXP nu) {
    if (LENGTH(mean) != LENGTH(nu))
      error("'mean' and 'nu' must be vectors of same length.");
    SEXP ans = PROTECT(allocVector(REALSXP, LENGTH(mean)));
    for(int i=0; i<LENGTH(mean); i++)
      REAL(ans)[i] = glmmtmb::compois_calc_var(REAL(mean)[i], REAL(nu)[i]);
    UNPROTECT(1);
    return ans;
  }
}

/* Quantile functions needed to simulate from truncated distributions */
extern "C" {
  double Rf_qnbinom(double p, double size, double prob, int lower_tail, int log_p);
  double Rf_qpois(double p, double lambda, int lower_tail, int log_p);
}

enum valid_family {
  gaussian_family = 0,
  binomial_family = 100,
  betabinomial_family =101,
  beta_family =200,
  Gamma_family =300,
  poisson_family =400,
  truncated_poisson_family =401,
  genpois_family =402,
  compois_family =403,
  truncated_genpois_family =404,
  truncated_compois_family =405,
  nbinom1_family =500,
  nbinom2_family =501,
  truncated_nbinom1_family =502,
  truncated_nbinom2_family =503,
  t_family =600,
  tweedie_family = 700
};

enum valid_link {
  log_link                 = 0,
  logit_link               = 1,
  probit_link              = 2,
  inverse_link             = 3,
  cloglog_link             = 4,
  identity_link            = 5
};

enum valid_covStruct {
  diag_covstruct = 0,
  us_covstruct   = 1,
  cs_covstruct   = 2,
  ar1_covstruct  = 3,
  ou_covstruct   = 4,
  exp_covstruct = 5,
  gau_covstruct = 6,
  mat_covstruct = 7,
  toep_covstruct = 8
};

enum valid_ziPredictCode {
  corrected_zipredictcode = 0,
  uncorrected_zipredictcode = 1,
  prob_zipredictcode = 2
};

template<class Type>
Type inverse_linkfun(Type eta, int link) {
  Type ans;
  switch (link) {
  case log_link:
    ans = exp(eta);
    break;
  case identity_link:
    ans = eta;
    break;
  case logit_link:
    ans = invlogit(eta);
    break;
  case probit_link:
    ans = pnorm(eta);
    break;
  case cloglog_link:
    ans = Type(1) - exp(-exp(eta));
    break;
  case inverse_link:
    ans = Type(1) / eta;
    break;
    // TODO: Implement remaining links
  default:
    error("Link not implemented!");
  } // End switch
  return ans;
}

/* logit transformed inverse_linkfun without losing too much
   accuracy */
template<class Type>
Type logit_inverse_linkfun(Type eta, int link) {
  Type ans;
  switch (link) {
  case logit_link:
    ans = eta;
    break;
  case probit_link:
    ans = glmmtmb::logit_pnorm(eta);
    break;
  case cloglog_link:
    ans = glmmtmb::logit_invcloglog(eta);
    break;
  default:
    ans = logit( inverse_linkfun(eta, link) );
  } // End switch
  return ans;
}

/* log transformed inverse_linkfun without losing too much accuracy */
template<class Type>
Type log_inverse_linkfun(Type eta, int link) {
  Type ans;
  switch (link) {
  case log_link:
    ans = eta;
    break;
  default:
    ans = log( inverse_linkfun(eta, link) );
  } // End switch
  return ans;
}

template <class Type>
struct per_term_info {
  // Input from R
  int blockCode;     // Code that defines structure
  int blockSize;     // Size of one block
  int blockReps;     // Repeat block number of times
  int blockNumTheta; // Parameter count per block
  matrix<Type> dist;
  vector<Type> times;// For ar1 case
  // Report output
  matrix<Type> corr;
  vector<Type> sd;
};

template <class Type>
struct terms_t : vector<per_term_info<Type> > {
  terms_t(SEXP x){
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP y = VECTOR_ELT(x, i);    // y = x[[i]]
      int blockCode = (int) REAL(getListElement(y, "blockCode", &isNumericScalar))[0];
      int blockSize = (int) REAL(getListElement(y, "blockSize", &isNumericScalar))[0];
      int blockReps = (int) REAL(getListElement(y, "blockReps", &isNumericScalar))[0];
      int blockNumTheta = (int) REAL(getListElement(y, "blockNumTheta", &isNumericScalar))[0];
      (*this)(i).blockCode = blockCode;
      (*this)(i).blockSize = blockSize;
      (*this)(i).blockReps = blockReps;
      (*this)(i).blockNumTheta = blockNumTheta;
      // Optionally, pass time vector:
      SEXP t = getListElement(y, "times");
      if(!isNull(t)){
	RObjectTestExpectedType(t, &isNumeric, "times");
	(*this)(i).times = asVector<Type>(t);
      }
      // Optionally, pass distance matrix:
      SEXP d = getListElement(y, "dist");
      if(!isNull(d)){
	RObjectTestExpectedType(d, &isMatrix, "dist");
	(*this)(i).dist = asMatrix<Type>(d);
      }
    }
  }
};

template <class Type>
Type termwise_nll(array<Type> &U, vector<Type> theta, per_term_info<Type>& term, bool do_simulate = false) {
  Type ans = 0;
  if (term.blockCode == diag_covstruct){
    // case: diag_covstruct
    vector<Type> sd = exp(theta);
    for(int i = 0; i < term.blockReps; i++){
      ans -= dnorm(vector<Type>(U.col(i)), Type(0), sd, true).sum();
      if (do_simulate) {
        U.col(i) = rnorm(Type(0), sd);
      }
    }
    term.sd = sd; // For report
  }
  else if (term.blockCode == us_covstruct){
    // case: us_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    vector<Type> corr_transf = theta.tail(theta.size() - n);
    vector<Type> sd = exp(logsd);
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
    density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  else if (term.blockCode == cs_covstruct){
    // case: cs_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    Type corr_transf = theta(n);
    vector<Type> sd = exp(logsd);
    Type a = Type(1) / (Type(n) - Type(1));
    Type rho = invlogit(corr_transf) * (Type(1) + a) - a;
    matrix<Type> corr(n,n);
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
	corr(i,j) = (i==j ? Type(1) : rho);
    density::MVNORM_t<Type> nldens(corr);
    density::VECSCALE_t<density::MVNORM_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  else if (term.blockCode == toep_covstruct){
    // case: toep_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    vector<Type> sd = exp(logsd);
    vector<Type> parms = theta.tail(n-1);              // Corr parms
    parms = parms / sqrt(Type(1.0) + parms * parms );  // Now in (-1,1)
    matrix<Type> corr(n,n);
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        corr(i,j) = (i==j ? Type(1) :
                     parms( (i > j ? i-j : j-i) - 1 ) );
    density::MVNORM_t<Type> nldens(corr);
    density::VECSCALE_t<density::MVNORM_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  else if (term.blockCode == ar1_covstruct){
    // case: ar1_covstruct
    //  * NOTE: Valid parameter space is phi in [-1, 1]
    //  * NOTE: 'times' not used as we assume unit distance between consecutive time points.
    int n = term.blockSize;
    Type logsd = theta(0);
    Type corr_transf = theta(1);
    Type phi = corr_transf / sqrt(1.0 + pow(corr_transf, 2));
    Type sd = exp(logsd);
    for(int j = 0; j < term.blockReps; j++){
      ans -= dnorm(U(0, j), Type(0), sd, true);   // Initialize
      if (do_simulate) {
        U(0, j) = rnorm(Type(0), sd);
      }
      for(int i=1; i<n; i++){
	ans -= dnorm(U(i, j), phi * U(i-1, j), sd * sqrt(1 - phi*phi), true);
        if (do_simulate) {
          U(i, j) = rnorm( phi * U(i-1, j), sd * sqrt(1 - phi*phi) );
        }
      }
    }
    // For consistency with output for other structs we report entire
    // covariance matrix.
    if(isDouble<Type>::value) { // Disable AD for this part
      term.corr.resize(n,n);
      term.sd.resize(n);
      for(int i=0; i<n; i++){
	term.sd(i) = sd;
	for(int j=0; j<n; j++){
	  term.corr(i,j) = pow(phi, abs(i-j));
	}
      }
    }
  }
  else if (term.blockCode == ou_covstruct){
    // case: ou_covstruct
    //  * NOTE: this is the continuous time version of ar1.
    //          One-step correlation must be non-negative
    //  * NOTE: 'times' assumed sorted !
    int n = term.times.size();
    Type logsd = theta(0);
    Type corr_transf = theta(1);
    Type sd = exp(logsd);
    for(int j = 0; j < term.blockReps; j++){
      ans -= dnorm(U(0, j), Type(0), sd, true);   // Initialize
      if (do_simulate) {
        U(0, j) = rnorm(Type(0), sd);
      }
      for(int i=1; i<n; i++){
	Type rho = exp(-exp(corr_transf) * (term.times(i) - term.times(i-1)));
	ans -= dnorm(U(i, j), rho * U(i-1, j), sd * sqrt(1 - rho*rho), true);
        if (do_simulate) {
          U(i, j) = rnorm( rho * U(i-1, j), sd * sqrt(1 - rho*rho));
        }
      }
    }
    // For consistency with output for other structs we report entire
    // covariance matrix.
    if(isDouble<Type>::value) { // Disable AD for this part
      term.corr.resize(n,n);
      term.sd.resize(n);
      for(int i=0; i<n; i++){
	term.sd(i) = sd;
	for(int j=0; j<n; j++){
	  term.corr(i,j) =
	    exp(-exp(corr_transf) * CppAD::abs(term.times(i) - term.times(j)));
	}
      }
    }
  }
  // Spatial correlation structures
  else if (term.blockCode == exp_covstruct ||
           term.blockCode == gau_covstruct ||
           term.blockCode == mat_covstruct){
    int n = term.blockSize;
    matrix<Type> dist = term.dist;
    if(! ( dist.cols() == n && dist.rows() == n ) )
      error ("Dimension of distance matrix must equal blocksize.");
    // First parameter is sd
    Type sd = exp( theta(0) );
    // Setup correlation matrix
    matrix<Type> corr(n,n);
    for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
        switch (term.blockCode) {
        case exp_covstruct:
          corr(i,j) = (i==j ? Type(1) : exp( -dist(i,j) * exp(-theta(1)) ) );
          break;
        case gau_covstruct:
          corr(i,j) = (i==j ? Type(1) : exp( -pow(dist(i,j),2) * exp(-2. * theta(1)) ) );
          break;
        case mat_covstruct:
          corr(i,j) = (i==j ? Type(1) : matern( dist(i,j),
                                                exp(theta(1)) /* range */,
                                                exp(theta(2)) /* smoothness */) );
          break;
        default:
          error("Not implemented");
        }
      }
    }
    density::MVNORM_t<Type> nldens(corr);
    density::SCALE_t<density::MVNORM_t<Type> > scnldens = density::SCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = corr;   // For report
    term.sd.resize(n);  // For report
    term.sd.fill(sd);
  }
  else error("covStruct not implemented!");
  return ans;
}

template <class Type>
Type allterms_nll(vector<Type> &u, vector<Type> theta,
		  vector<per_term_info<Type> >& terms,
                  bool do_simulate = false) {
  Type ans = 0;
  int upointer = 0;
  int tpointer = 0;
  int nr, np = 0, offset;
  for(int i=0; i < terms.size(); i++){
    nr = terms(i).blockSize * terms(i).blockReps;
    // Note: 'blockNumTheta=0' ==> Same parameters as previous term.
    bool emptyTheta = ( terms(i).blockNumTheta == 0 );
    offset = ( emptyTheta ? -np : 0 );
    np     = ( emptyTheta ?  np : terms(i).blockNumTheta );
    vector<int> dim(2);
    dim << terms(i).blockSize, terms(i).blockReps;
    array<Type> useg( &u(upointer), dim);
    vector<Type> tseg = theta.segment(tpointer + offset, np);
    ans += termwise_nll(useg, tseg, terms(i), do_simulate);
    upointer += nr;
    tpointer += terms(i).blockNumTheta;
  }
  return ans;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(X);
  DATA_SPARSE_MATRIX(Z);
  DATA_MATRIX(Xzi);
  DATA_SPARSE_MATRIX(Zzi);
  DATA_MATRIX(Xd);
  DATA_VECTOR(yobs);
  DATA_VECTOR(size); //only used in binomial
  DATA_VECTOR(weights);
  DATA_VECTOR(offset);
  DATA_VECTOR(zioffset);
  DATA_VECTOR(doffset);

  // Define covariance structure for the conditional model
  DATA_STRUCT(terms, terms_t);

  // Define covariance structure for the zero inflation
  DATA_STRUCT(termszi, terms_t);

  // Parameters related to design matrices
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(betazi);
  PARAMETER_VECTOR(b);
  PARAMETER_VECTOR(bzi);
  PARAMETER_VECTOR(betad);

  // Joint vector of covariance parameters
  PARAMETER_VECTOR(theta);
  PARAMETER_VECTOR(thetazi);

  // Extra family specific parameters (e.g. tweedie)
  PARAMETER_VECTOR(thetaf);

  DATA_INTEGER(family);
  DATA_INTEGER(link);

  // Flags
  DATA_INTEGER(ziPredictCode);
  bool zi_flag = (betazi.size() > 0);
  DATA_INTEGER(doPredict);
  DATA_IVECTOR(whichPredict);

  // One-Step-Ahead (OSA) residuals
  DATA_VECTOR_INDICATOR(keep, yobs);

  // Joint negative log-likelihood
  Type jnll = 0;

  // Random effects
  jnll += allterms_nll(b, theta, terms, this->do_simulate);
  jnll += allterms_nll(bzi, thetazi, termszi, this->do_simulate);

  // Linear predictor
  vector<Type> eta = X * beta + Z * b + offset;
  vector<Type> etazi = Xzi * betazi + Zzi * bzi + zioffset;
  vector<Type> etad = Xd * betad + doffset;

  // Apply link
  vector<Type> mu(eta.size());
  for (int i = 0; i < mu.size(); i++)
    mu(i) = inverse_linkfun(eta(i), link);
  vector<Type> pz = invlogit(etazi);
  vector<Type> phi = exp(etad);

  // Observation likelihood
  Type s1, s2, s3, log_nzprob;
  Type tmp_loglik;
  for (int i=0; i < yobs.size(); i++){
    if ( !glmmtmb::isNA(yobs(i)) ) {
      switch (family) {
      case gaussian_family:
        tmp_loglik = dnorm(yobs(i), mu(i), sqrt(phi(i)), true);
        SIMULATE{yobs(i) = rnorm(mu(i), sqrt(phi(i)));}
        break;
      case poisson_family:
        tmp_loglik = dpois(yobs(i), mu(i), true);
        SIMULATE{yobs(i) = rpois(mu(i));}
        break;
      case binomial_family:
        s1 = logit_inverse_linkfun(eta(i), link); // logit(p)
        tmp_loglik = dbinom_robust(yobs(i), size(i), s1, true);
        SIMULATE{yobs(i) = rbinom(size(i), mu(i));}
        break;
      case Gamma_family:
        s1 = phi(i);           // shape
        s2 = mu(i) / phi(i);   // scale
        tmp_loglik = dgamma(yobs(i), s1, s2, true);
        SIMULATE{yobs(i) = rgamma(s1, s2);}
        break;
      case beta_family:
        // parameterization after Ferrari and Cribari-Neto 2004, betareg package
        s1 = mu(i)*phi(i);
        s2 = (Type(1)-mu(i))*phi(i);
        tmp_loglik = dbeta(yobs(i), s1, s2, true);
        SIMULATE{yobs(i) = rbeta(s1, s2);}
        break;
      case betabinomial_family:
        s1 = mu(i)*phi(i); // s1 = mu(i) * mu(i) / phi(i);
        s2 = (Type(1)-mu(i))*phi(i); // phi(i) / mu(i);
        tmp_loglik = glmmtmb::dbetabinom(yobs(i), s1, s2, size(i), true);
        SIMULATE {
          yobs(i) = rbinom(size(i), rbeta(s1, s2) );
        }
        break;
      case nbinom1_family:
      case truncated_nbinom1_family:
        // Was:
        //   s1 = mu(i);
        //   s2 = mu(i) * (Type(1)+phi(i));  // (1+phi) guarantees that var >= mu
        //   tmp_loglik = dnbinom2(yobs(i), s1, s2, true);
        s1 = log_inverse_linkfun(eta(i), link);          // log(mu)
        s2 = s1 + etad(i) ;                              // log(var - mu)
        tmp_loglik = dnbinom_robust(yobs(i), s1, s2, true);
        SIMULATE {
          s1 = mu(i);
          s2 = mu(i) * (Type(1)+phi(i));  // (1+phi) guarantees that var >= mu
          yobs(i) = rnbinom2(s1, s2);
        }
        if( family == truncated_nbinom1_family ) {
          // s3 := log( 1. + phi(i) )
          s3 = logspace_add( Type(0), etad(i) );
          log_nzprob = logspace_sub( Type(0), -mu(i) / phi(i) * s3 ); // 1-prob(0)
          tmp_loglik -= log_nzprob;
          if( yobs(i) < Type(1) ) tmp_loglik = -INFINITY;
          SIMULATE{
            s1 = mu(i)/phi(i);//sz
            s2 = 1/(1+phi(i)); //pb
            yobs(i) = Rf_qnbinom(asDouble(runif(dnbinom(Type(0), s1, s2), Type(1))), asDouble(s1), asDouble(s2), 1, 0);
          }
        }
        break;
      case nbinom2_family:
      case truncated_nbinom2_family:
        // Was:
        //   s1 = mu(i);
        //   s2 = mu(i) * (Type(1) + mu(i) / phi(i));
        //   tmp_loglik = dnbinom2(yobs(i), s1, s2, true);
        s1 = log_inverse_linkfun(eta(i), link);          // log(mu)
        s2 = 2. * s1 - etad(i) ;                         // log(var - mu)
        tmp_loglik = dnbinom_robust(yobs(i), s1, s2, true);
        SIMULATE {
          s1 = mu(i);
          s2 = mu(i) * (Type(1) + mu(i) / phi(i));
          yobs(i) = rnbinom2(s1, s2);
        }
        if (family == truncated_nbinom2_family) {
          // s3 := log( 1. + mu(i) / phi(i) )
          s3         = logspace_add( Type(0), s1 - etad(i) );
          log_nzprob = logspace_sub( Type(0), -phi(i) * s3 );
          tmp_loglik -= log_nzprob;
          if( yobs(i) < Type(1) ) tmp_loglik = -INFINITY;
          SIMULATE{
            s1 = phi(i); //sz
            s2 = phi(i)/(phi(i)+mu(i)); //pb
            yobs(i) = Rf_qnbinom(asDouble(runif(dnbinom(Type(0), s1, s2), Type(1))), asDouble(s1), asDouble(s2), 1, 0);
          }
        }
        break;
      case truncated_poisson_family:
        // Was:
        //   if (mu(i)<1e-6) {
        //     nzprob = mu(i)*(1-mu(i)/2);
        //   } else {
        //     nzprob = 1-exp(-mu(i));
        //   }
        // log(nzprob) = log( 1 - exp(-mu(i)) )
        log_nzprob = logspace_sub(Type(0), -mu(i));
        tmp_loglik = dpois(yobs(i), mu(i), true) - log_nzprob;
        if( yobs(i) < Type(1) ) tmp_loglik = -INFINITY;
        SIMULATE{
          yobs(i) = Rf_qpois(asDouble(runif(dpois(Type(0), mu(i)), Type(1))), asDouble(mu(i)), 1, 0);
        }
        break;
     case genpois_family:
        s1 = mu(i) / sqrt(phi(i)); //theta
        s2 = Type(1) - Type(1)/sqrt(phi(i)); //lambda
        tmp_loglik = glmmtmb::dgenpois(yobs(i), s1, s2, true);
        SIMULATE{yobs(i)=glmmtmb::rgenpois(mu(i) / sqrt(phi(i)), Type(1) - Type(1)/sqrt(phi(i)));}
        break;
      case truncated_genpois_family:
        s1 = mu(i) / sqrt(phi(i)); //theta
        s2 = Type(1) - Type(1)/sqrt(phi(i)); //lambda
        log_nzprob = logspace_sub(Type(0), -s1);
        tmp_loglik = glmmtmb::dgenpois(yobs(i), s1, s2, true) - log_nzprob;
        if( yobs(i) < Type(1) ) tmp_loglik = -INFINITY;
        SIMULATE{yobs(i)=glmmtmb::rtruncated_genpois(mu(i) / sqrt(phi(i)), Type(1) - Type(1)/sqrt(phi(i)));}
        break;
      case compois_family:
        s1 = mu(i); //mean
        s2 = 1/phi(i); //nu
        tmp_loglik = dcompois2(yobs(i), s1, s2, true);
        SIMULATE{yobs(i)=rcompois2(mu(i), 1/phi(i));}
        break;
      case truncated_compois_family:
        s1 = mu(i); //mean
        s2 = 1/phi(i); //nu
        log_nzprob = logspace_sub(Type(0), dcompois2(Type(0), s1, s2, true));
        tmp_loglik = dcompois2(yobs(i), s1, s2, true) - log_nzprob;
        if( yobs(i) < Type(1) ) tmp_loglik = -INFINITY;
        SIMULATE{yobs(i)=glmmtmb::rtruncated_compois2(mu(i), 1/phi(i));}
        break;
      case tweedie_family:
        s1 = mu(i);  // mean
        s2 = phi(i); // phi
        s3 = invlogit(thetaf(0)) + Type(1); // p, 1<p<2
        tmp_loglik = dtweedie(yobs(i), s1, s2, s3, true);
        SIMULATE {
          yobs(i) = glmmtmb::rtweedie(s1, s2, s3);
        }
        break;
      default:
        error("Family not implemented!");
      } // End switch

      // Add zero inflation
      if(zi_flag){
        Type logit_pz = etazi(i) ;
        Type log_pz   = -logspace_add( Type(0) , -logit_pz );
        Type log_1mpz = -logspace_add( Type(0) ,  logit_pz );
        if(yobs(i) == Type(0)){
          // Was:
          //   tmp_loglik = log( pz(i) + (1.0 - pz(i)) * exp(tmp_loglik) );
          tmp_loglik = logspace_add( log_pz, log_1mpz + tmp_loglik );
        } else {
          // Was:
          //   tmp_loglik += log( 1.0 - pz(i) );
          tmp_loglik += log_1mpz ;
        }
        SIMULATE{yobs(i) = yobs(i)*rbinom(Type(1), Type(1)-pz(i));}
      }
      tmp_loglik *= weights(i);
      // Add up
      jnll -= keep(i) * tmp_loglik;
    }
  }

  // Report / ADreport / Simulate Report
  vector<matrix<Type> > corr(terms.size());
  vector<vector<Type> > sd(terms.size());
  for(int i=0; i<terms.size(); i++){
    // NOTE: Dummy terms reported as empty
    if(terms(i).blockNumTheta > 0){
      corr(i) = terms(i).corr;
      sd(i) = terms(i).sd;
    }
  }
  vector<matrix<Type> > corrzi(termszi.size());
  vector<vector<Type> > sdzi(termszi.size());
  for(int i=0; i<termszi.size(); i++){
    // NOTE: Dummy terms reported as empty
    if(termszi(i).blockNumTheta > 0){
      corrzi(i) = termszi(i).corr;
      sdzi(i) = termszi(i).sd;
    }
  }

  REPORT(corr);
  REPORT(sd);
  REPORT(corrzi);
  REPORT(sdzi);
  SIMULATE {
    REPORT(yobs);
    REPORT(b);
    REPORT(bzi);
  }
  // For predict
  if(zi_flag) {
    switch(ziPredictCode){
    case corrected_zipredictcode:
      mu *= (Type(1) - pz); // Account for zi in prediction
      break;
    case uncorrected_zipredictcode:
      //mu = mu; // Predict mean of 'family' //comented out for clang 7.0.0. with no effect
      break;
    case prob_zipredictcode:
      mu = pz; // Predict zi probability
      break;
    default:
      error("Invalid 'ziPredictCode'");
    }
  }

  whichPredict -= 1; // R-index -> C-index
  vector<Type> mu_predict = mu(whichPredict);
  REPORT(mu_predict);
  // ADREPORT expensive for long vectors - only needed by predict()
  // method.
  if (doPredict) ADREPORT(mu_predict);

  return jnll;
}

