#include <TMB.hpp>
#include "init.h"
#include "distrib.h"

// don't need to include omp.h; we get it via TMB.hpp

namespace glmmtmb{
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
bool notFinite(Type x) {
  return (!R_FINITE(asDouble(x)));
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
  identity_link            = 5,
  sqrt_link                = 6
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
  toep_covstruct = 8,
  rr_covstruct = 9
};

enum valid_ziPredictCode {
  corrected_zipredictcode = 0,
  uncorrected_zipredictcode = 1,
  prob_zipredictcode = 2,
  disp_zipredictcode = 3
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
  case sqrt_link:
    ans = eta*eta; // pow(eta, Type(2)) doesn't work ... ?
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
  case logit_link:
    ans = -logspace_add(Type(0), -eta);
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
  matrix<Type> fact_load; // For rr case
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
  else if (term.blockCode == rr_covstruct){
    // case: reduced rank
    for(int i = 0; i < term.blockReps; i++){
      ans -= dnorm(vector<Type>(U.col(i)), Type(0), 1, true).sum();
      if (do_simulate) {
        U.col(i) = rnorm(U.rows(), Type(0), Type(1));
      }
    }

    int p = term.blockSize;
    int nt = theta.size();
    int rank = (2*p + 1 -  (int)sqrt(pow(2.0*p + 1, 2) - 8*nt) ) / 2 ;
    matrix<Type> Lambda(p, rank);
    vector<Type> lam_diag = theta.head(rank);
    vector<Type> lam_lower = theta.tail(nt - rank);
    for (int j = 0; j < rank; j++){
      for (int i = 0; i < p; i++){
        if (j > i)
          Lambda(i, j) = 0;
        else if(i == j)
          Lambda(i, j) = lam_diag(j);
        else
          Lambda(i, j) = lam_lower(j*p - (j + 1)*j/2 + i - 1 - j); //Fills by column
      }
    }

    for(int i = 0; i < term.blockReps; i++){
      vector<Type> usub = U.col(i).segment(0, rank);
      U.col(i) = Lambda * usub;
    }

    term.fact_load = Lambda;
    if(isDouble<Type>::value) {
      term.corr = Lambda * Lambda.transpose();
      term.sd = term.corr.diagonal().array().sqrt();
      term.corr.array() /= term.sd * term.sd.transpose();
    }
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
#ifdef _OPENMP
  this -> max_parallel_regions = omp_get_max_threads();
  // std::cout << "OpenMP max_parallel_regions=" << this -> max_parallel_regions << "\n";
#else
  this -> max_parallel_regions = 1;
  // std::cout << "no OpenMP (max_parallel_regions=1)\n";
#endif

  DATA_MATRIX(X);
  bool sparseX = X.rows()==0 && X.cols()==0;
  DATA_SPARSE_MATRIX(Z);
  DATA_MATRIX(Xzi);
  bool sparseXzi = Xzi.rows()==0 && Xzi.cols()==0;
  DATA_SPARSE_MATRIX(Zzi);
  DATA_MATRIX(Xd);
  bool sparseXd = Xd.rows()==0 && Xd.cols()==0;
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
  Type jnll=0;

  // Random effects
  PARALLEL_REGION jnll += allterms_nll(b, theta, terms, this->do_simulate);
  PARALLEL_REGION jnll += allterms_nll(bzi, thetazi, termszi, this->do_simulate);

  // Linear predictor
  vector<Type> eta = Z * b + offset;
  if (!sparseX) {
    eta += X*beta;
  } else {
    DATA_SPARSE_MATRIX(XS);
    eta += XS*beta;
  }
  vector<Type> etazi = Zzi * bzi + zioffset;
  if (!sparseXzi) {
    etazi += Xzi*betazi;
  } else {
    DATA_SPARSE_MATRIX(XziS);
    etazi += XziS*betazi;
  }
  vector<Type> etad = doffset;
  if (!sparseXd) {
    etad += Xd*betad;
  } else {
    DATA_SPARSE_MATRIX(XdS);
    etad += XdS*betad;
  }

  // Apply link
  vector<Type> mu(eta.size());
  for (int i = 0; i < mu.size(); i++)
    mu(i) = inverse_linkfun(eta(i), link);
  vector<Type> pz = invlogit(etazi);
  vector<Type> phi = exp(etad);

// "zero-truncated" likelihood: ignore zeros in positive distributions
// exact zero: use for positive distributions (Gamma, beta)
#define zt_lik_zero(x,loglik_exp) (zi_flag && (x == Type(0)) ? -INFINITY : loglik_exp)
// close to zero: use for count data (cf binomial()$initialize)
#define zt_lik_nearzero(x,loglik_exp) (zi_flag && (x < Type(0.001)) ? -INFINITY : loglik_exp)

  // Observation likelihood
  Type s1, s2, s3, log_nzprob;
  Type tmp_loglik;
  for (int i=0; i < yobs.size(); i++) PARALLEL_REGION {
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
        tmp_loglik = zt_lik_zero(yobs(i),dgamma(yobs(i), s1, s2, true));
        SIMULATE{yobs(i) = rgamma(s1, s2);}
        break;
      case beta_family:
        // parameterization after Ferrari and Cribari-Neto 2004, betareg package
        s1 = mu(i)*phi(i);
        s2 = (Type(1)-mu(i))*phi(i);
        tmp_loglik = zt_lik_zero(yobs(i),dbeta(yobs(i), s1, s2, true));
        SIMULATE{yobs(i) = rbeta(s1, s2);}
        break;
      case betabinomial_family:
        // Transform to logit scale independent of link
        s3 = logit_inverse_linkfun(eta(i), link); // logit(p)
        // Was: s1 = mu(i) * phi(i);
        s1 = log_inverse_linkfun( s3, logit_link) + log(phi(i)); // s1 = log(mu*phi)
        // Was: s2 = (Type(1) - mu(i)) * phi(i);
        s2 = log_inverse_linkfun(-s3, logit_link) + log(phi(i)); // s2 = log((1-mu)*phi)
        tmp_loglik = glmmtmb::dbetabinom_robust(yobs(i), s1, s2, size(i), true);
        SIMULATE {
          yobs(i) = rbinom(size(i), rbeta(exp(s1), exp(s2)) );
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
	if (family != truncated_nbinom1_family) {
		SIMULATE {
			s1 = mu(i);
			s2 = mu(i) * (Type(1)+phi(i));  // (1+phi) guarantees that var >= mu
			yobs(i) = rnbinom2(s1, s2);
		}
	} else {
          s3 = logspace_add( Type(0), etad(i) );                // log(1. + phi(i)
          log_nzprob = logspace_sub( Type(0), -mu(i) / phi(i) * s3 ); // 1-prob(0)
          tmp_loglik -= log_nzprob;
	  tmp_loglik = zt_lik_nearzero(yobs(i), tmp_loglik);
          SIMULATE{
            s1 = mu(i)/phi(i); //sz
	    yobs(i) = glmmtmb::rtruncated_nbinom(asDouble(s1), 0, asDouble(mu(i)));
          }
        }
        break;
      case nbinom2_family:
      case truncated_nbinom2_family:
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
          tmp_loglik = zt_lik_nearzero( yobs(i), tmp_loglik);
          SIMULATE{
		  yobs(i) = glmmtmb::rtruncated_nbinom(asDouble(phi(i)), 0, asDouble(mu(i)));
          }
        }
        break;
      case truncated_poisson_family:
        log_nzprob = logspace_sub(Type(0), -mu(i));  // log(1-exp(-mu(i))) = P(x>0)
        tmp_loglik = dpois(yobs(i), mu(i), true) - log_nzprob;
        tmp_loglik = zt_lik_nearzero(yobs(i), tmp_loglik);
        SIMULATE{
		yobs(i) = glmmtmb::rtruncated_poisson(0, asDouble(mu(i)));
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
        tmp_loglik = zt_lik_nearzero(yobs(i),
		    glmmtmb::dgenpois(yobs(i), s1, s2, true) - log_nzprob);
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
        tmp_loglik = zt_lik_nearzero(yobs(i),
			    dcompois2(yobs(i), s1, s2, true) - log_nzprob);
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

  vector<matrix<Type> > fact_load(terms.size());
  for(int i=0; i<terms.size(); i++){
    // NOTE: Dummy terms reported as empty
    if(terms(i).blockNumTheta > 0){
      fact_load(i) = terms(i).fact_load;
    }
  }

  REPORT(corr);
  REPORT(sd);
  REPORT(corrzi);
  REPORT(sdzi);
  REPORT(fact_load);
  SIMULATE {
    REPORT(yobs);
    REPORT(b);
    REPORT(bzi);
  }
  // For predict
  if(ziPredictCode==disp_zipredictcode){ //zi irrelevant; just reusing variable
    switch(family){
    case gaussian_family:
      mu = sqrt(phi);
      break;
    case Gamma_family:
      mu = 1/sqrt(phi);
      break;
    default:
      mu = phi;
    }
  } else {
  if(zi_flag) {
    switch(ziPredictCode){
    case corrected_zipredictcode:
      mu *= (Type(1) - pz); // Account for zi in prediction
      break;
    case uncorrected_zipredictcode:
      //mu = mu; // Predict mean of 'family' //comented out for clang 7.0.0. with no effect
      break;
    case prob_zipredictcode:
      mu = pz;     // Predict zi probability
      eta = etazi; // want to return linear pred for zi
      break;
    default:
      error("Invalid 'ziPredictCode'");
    }
  }}

  whichPredict -= 1; // R-index -> C-index
  vector<Type> mu_predict = mu(whichPredict);
  vector<Type> eta_predict = eta(whichPredict);
  REPORT(mu_predict);
  REPORT(eta_predict);
  // ADREPORT expensive for long vectors - only needed by predict() method
  if (doPredict==1) {
	  ADREPORT(mu_predict);
  } else if (doPredict == 2) {
	  ADREPORT(eta_predict);
  }

  return jnll;
}



