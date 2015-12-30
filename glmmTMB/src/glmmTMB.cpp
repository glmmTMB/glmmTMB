#include <TMB.hpp>

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
  bool isNA(Type x){
    return R_IsNA(asDouble(x));
  }
}

enum valid_family {
  gaussian_family = 0,
  binomial_family = 100,
  betabinomial_family =101,
  beta_family =200,
  Gamma_family =300,
  poisson_family =400,
  truncated_poisson_family =401,
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
  ar1_covstruct  = 3
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
    }
  }
};

template <class Type>
Type termwise_nll(vector<Type> u, vector<Type> theta, per_term_info<Type>& term) {
  Type ans = 0;
  vector<int> dim(2);
  dim << term.blockSize, term.blockReps;
  array<Type> U(u, dim); // Note: Fill columnwise
  
  if (term.blockCode == diag_covstruct){
    // case: diag_covstruct
    vector<Type> sd = exp(theta);
    for(int i = 0; i < term.blockReps; i++){
      ans -= dnorm(vector<Type>(U.col(i)), Type(0), sd, true).sum();
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
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  else if (term.blockCode == ar1_covstruct){
    // case: ar1_covstruct
    //  * NOTE: 'times' assumed sorted !
    int n = term.times.size();
    Type logsd = theta(0);
    Type corr_transf = theta(1);
    Type sd = exp(logsd);
    for(int j = 0; j < term.blockReps; j++){
      ans -= dnorm(U(0, j), Type(0), sd, true);   // Initialize
      for(int i=1; i<n; i++){
	Type rho = exp(-exp(corr_transf) * (term.times(i) - term.times(i-1)));
	ans -= dnorm(U(i, j), rho * U(i-1, j), sd * sqrt(1 - rho*rho), true);
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
	    exp(-exp(corr_transf) * abs(term.times(i) - term.times(j)));
	}
      }
    }
  }
  else error("covStruct not implemented!");
  return ans;
}

template <class Type>
Type allterms_nll(vector<Type> u, vector<Type> theta, 
		  vector<per_term_info<Type> >& terms) {
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
    vector<Type> useg =     u.segment(upointer         , nr);
    vector<Type> tseg = theta.segment(tpointer + offset, np);
    ans += termwise_nll(useg, tseg, terms(i));
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
  DATA_VECTOR(weights);
  DATA_VECTOR(offset);

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

  DATA_INTEGER(family);
  DATA_INTEGER(link);

  // Flags
  DATA_INTEGER(ziPredictCode);
  bool zi_flag = (betazi.size() > 0);
  DATA_INTEGER(doPredict);

  // Joint negative log-likelihood
  Type jnll = 0;

  // Random effects
  jnll += allterms_nll(b, theta, terms);
  jnll += allterms_nll(bzi, thetazi, termszi);

  // Linear predictor
  vector<Type> eta = X * beta + Z * b + offset;
  vector<Type> etazi = Xzi * betazi + Zzi * bzi;
  vector<Type> etad = Xd * betad;

  // Apply link
  vector<Type> mu(eta.size());
  for (int i = 0; i < mu.size(); i++)
    mu(i) = inverse_linkfun(eta(i), link);
  vector<Type> pz = invlogit(etazi);
  vector<Type> phi = exp(etad);

  // Observation likelihood
  Type s1, s2, stmp;
  Type tmp_loglik;
  for (int i=0; i < yobs.size(); i++){
    if ( !glmmtmb::isNA(yobs(i)) ) {
      switch (family) {
      case gaussian_family:
	tmp_loglik = weights(i) * dnorm(yobs(i), mu(i), sqrt(phi(i)), true);
	break;
      case poisson_family:
	tmp_loglik = weights(i) * dpois(yobs(i), mu(i), true);
	break;
      case binomial_family:
	tmp_loglik = dbinom(yobs(i) * weights(i), weights(i), mu(i), true);
	break;
      case Gamma_family:
	s1 = phi(i);           // shape
	s2 = mu(i) / phi(i);   // scale
	tmp_loglik = weights(i) * dgamma(yobs(i), s1, s2, true);
	break;
      case beta_family:
        // parameterization after Ferrari and Cribari-Neto 2004, betareg package
        s1 = mu(i)*phi(i);
        s2 = (Type(1)-mu(i))*phi(i);
	tmp_loglik = weights(i) * dbeta(yobs(i), s1, s2, true);
	break;
      case betabinomial_family:
	s1 = mu(i) * mu(i) / phi(i);
	s2 = phi(i) / mu(i);
	tmp_loglik = glmmtmb::dbetabinom(yobs(i) * weights(i), s1, s2, weights(i), true);
	break;
      case nbinom1_family:
	s1 = mu(i);
	s2 = mu(i) * (Type(1)+phi(i));  // (1+phi) guarantees that var >= mu
	tmp_loglik = weights(i) * dnbinom2(yobs(i), s1, s2, true);
	break;
      case nbinom2_family:
	s1 = mu(i);
	s2 = mu(i) * (Type(1) + mu(i) / phi(i));
	tmp_loglik = weights(i) * dnbinom2(yobs(i), s1, s2, true);
	break;
	// TODO: Implement remaining families
      default:
	error("Family not implemented!");
      } // End switch

      // Add zero inflation
      if(zi_flag){
	if(yobs(i) == Type(0)){
	  tmp_loglik = log( pz(i) + (1.0 - pz(i)) * exp(tmp_loglik) );
	} else {
	  tmp_loglik += log( 1.0 - pz(i) );
	}
      }

      // Add up
      jnll -= tmp_loglik;
    }
  }

  // Report / ADreport
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

  // For predict
  if(zi_flag) {
    switch(ziPredictCode){
    case corrected_zipredictcode:
      mu *= (Type(1) - pz); // Account for zi in prediction
      break;
    case uncorrected_zipredictcode:
      mu = mu; // Predict mean of 'family'
      break;
    case prob_zipredictcode:
      mu = pz; // Predict zi probability
      break;
    default:
      error("Invalid 'ziPredictCode'");
    }
  }

  REPORT(mu);
  // ADREPORT expensive for long vectors - only needed by predict()
  // method. FIXME: May even consider reducing mu to some subset
  // before ADREPORTing.
  if (doPredict) ADREPORT(mu);

  return jnll;
}

