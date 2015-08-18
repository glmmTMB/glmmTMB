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
  diag_covstruct = 0,
  us_covstruct   = 1,
  cs_covstruct   = 2,
  ar1_covstruct  = 3
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
  vector<Type> sd;
  switch (blockCode){
  case diag_covstruct:
    sd = exp(theta);
    for(int i = 0; i < blockReps; i++){
      ans -= dnorm(vector<Type>(U.col(i)), Type(0), sd, true).sum();
    }
    break;
  case us_covstruct:
    break;
  default:
    error("covStruct not implemented!");
  } // End switch
  return ans;
}

template <class Type>
Type allterms_nll(vector<Type> u, vector<Type> theta, vector<int> blockCode, vector<int> blockSize, vector<int> blockReps, vector<int> blockNumTheta){
  Type ans = 0;
  int upointer = 0;
  int tpointer = 0;
  int nr;
  for(int i=0; i < blockCode.size(); i++){
    nr = blockSize(i) * blockReps(i);
    vector<Type> useg = u.segment(upointer, nr);
    vector<Type> tseg = theta.segment(tpointer, blockNumTheta(i));
    ans += termwise_nll(useg, tseg, blockCode(i), blockSize(i), blockReps(i), blockNumTheta(i));
    upointer += nr;
    tpointer += blockNumTheta(i);
  }
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
  PARAMETER_VECTOR(theta);
  PARAMETER_VECTOR(thetazi);

  DATA_INTEGER(family);
  DATA_INTEGER(link);

  // Flags
  bool zi_flag = (betazi.size() > 0);

  // Joint negative log-likelihood
  Type jnll = 0;

  // Random effects
  jnll += allterms_nll(b, theta, blockCode, blockSize, blockReps, blockNumTheta);
  jnll += allterms_nll(bzi, thetazi, blockCodezi, blockSizezi, blockRepszi, blockNumThetazi);

  // Linear predictor
  vector<Type> eta = X * beta + Z * b;
  vector<Type> etazi = Xzi * betazi + Zzi * bzi;
  vector<Type> etad = Xd * betad;

  // Apply link
  vector<Type> mu(eta.size());
  for (int i = 0; i < mu.size(); i++)
    mu(i) = inverse_linkfun(eta(i), link);
  vector<Type> pz = invlogit(etazi);
  vector<Type> phi = exp(etad);

  // Observation likelihood
  Type tmp_loglik;
  for (int i=0; i < yobs.size(); i++){

    switch (family) {
    case gaussian_family:
      tmp_loglik = dnorm(yobs(i), mu(i), sqrt(phi(i)), true);
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

