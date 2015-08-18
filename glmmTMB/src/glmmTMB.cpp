#include <TMB.hpp>

enum valid_family {
  gaussian_family = 0,
  binomial_family = 100,
  betabinomial_family =101,
  beta_family =200,
  gamma_family =300,
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
  diag = 0,
  us   = 1,
  cs   = 2,
  ar1  = 3
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
Type termwise_nll(vector<Type> u, vector<Type> theta, int blockCode, int blockSize, int blockReps, int blockNumTheta){
  Type ans = 0;
  vector<int> dim(2);
  dim << blockSize, blockReps;
  array<Type> U(u, dim); // Note: Fill columnwise
  switch (blockCode){
  case diag:
    vector<Type> sd = exp(theta);
    for(int i = 0; i < blockReps; i++){
      ans -= dnorm(vector<Type>(U.col(i)), Type(0), sd, true).sum();
    }
    break;
  case us:

  default:
    error("covStruct not implemented!");
  } // End switch
  return ans;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(X);
  DATA_MATRIX(Z);
  DATA_MATRIX(Xzi);
  DATA_MATRIX(Zzi);
  DATA_MATRIX(Xd);
  DATA_VECTOR(yobs);

  // Define covariance structure for the conditional model
  DATA_IVECTOR(blockCode);     // Code that defines structure
  DATA_IVECTOR(blockSize);     // Size of one block
  DATA_IVECTOR(blockReps);     // Repeat block number of times
  DATA_IVECTOR(blockNumTheta); // Parameter count per block

  // Define covariance structure for the zero inflation
  DATA_IVECTOR(blockCodezi);     // Code that defines structure
  DATA_IVECTOR(blockSizezi);     // Size of one block
  DATA_IVECTOR(blockRepszi);     // Repeat block number of times
  DATA_IVECTOR(blockNumThetazi); // Parameter count per block

  // Parameters related to design matrices
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(betazi); // TODO: Notation ?
  PARAMETER_VECTOR(b);
  PARAMETER_VECTOR(bzi);
  PARAMETER_VECTOR(betad);

  // Joint vector of covariance parameters
  PARAMETER(theta);
  PARAMETER(thetazi);

  DATA_INTEGER(family);
  DATA_INTEGER(link);

  // Joint negative log-likelihood
  Type jnll = 0;

  // Random effects
  // TODO: Other corr structures
  if (b.size() > 0) {
    // vector<Type> sigma = exp(logsigma);
    // jnll -= dnorm(b, Type(0), sigma, true).sum();
  }

  // Linear predictor
  vector<Type> eta = X * beta + Z * b;

  // Apply link
  vector<Type> mu(eta.size());
  for (int i = 0; i < mu.size(); i++)
    mu(i) = inverse_linkfun(eta(i), link);

  // Zero-inflation
  // TODO: Allow multiple links for zero-infl as well ?
  bool zi_flag = (betazi.size() > 0);
  vector<Type> pz;
  if(zi_flag) {
    vector<Type> etazi = Zzi * betazi;
    pz = invlogit(etazi);
  }

  // Observation likelihood
  Type tmp_loglik;
  for (int i=0; i<yobs.size(); i++){

    switch (family) {
    case gaussian_family:
      // tmp_loglik = dnorm(yobs(i), mu(i), exp(logalpha), true);
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

