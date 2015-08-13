#include <TMB.hpp>

enum valid_family {
  poisson_family           = 0,
  binomial_family          = 1,
  negative_binomial_family = 2,
  Gamma_family             = 3,
  beta_family              = 4,
  gaussian_family          = 5,
  truncated_poisson_family = 6,
  trunc_NB_family          = 7,
  logistic_family          = 8,
  betabinomial_family      = 9
};

enum valid_link {
  log_link                 = 0,
  logit_link               = 1,
  probit_link              = 2,
  inverse_link             = 3,
  cloglog_link             = 4,
  identity_link            = 5
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
    // TODO: Implement remaining links
  default:
    error("Link not implemented!");
  } // End switch
  return ans;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(X);
  DATA_MATRIX(Z);
  DATA_MATRIX(ZI); // TODO: Notation ?
  DATA_VECTOR(yobs); // TODO: Pass matrix ?

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(b);
  PARAMETER_VECTOR(betazi); // TODO: Notation ?
  PARAMETER_VECTOR(logsigma);
  PARAMETER(logalpha); // TODO: Notation (taken from glmmadmb)

  DATA_INTEGER(family);
  DATA_INTEGER(link);

  // Joint negative log-likelihood
  Type jnll = 0;

  // Random effects
  // TODO: Other corr structures
  if (b.size() > 0) {
    vector<Type> sigma = exp(logsigma);
    jnll -= dnorm(b, Type(0), sigma, true).sum();
  }

  // Linear predictor
  vector<Type> eta = X * beta + Z * b;

  // Apply link
  vector<Type> mu = inverse_linkfun(eta, link);

  // Zero-inflation
  // TODO: Allow multiple links for zero-infl as well ?
  bool zi_flag = (betazi.size() > 0);
  vector<Type> pz;
  if(zi_flag) {
    vector<Type> etazi = ZI * betazi;
    pz = invlogit(etazi);
  }

  // Observation likelihood
  Type tmp_loglik;
  for (int i=0; i<yobs.size(); i++){

    switch (family) {
    case gaussian_family:
      tmp_loglik = dnorm(yobs(i), mu(i), exp(logalpha), true);
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

  return jnll;
}

