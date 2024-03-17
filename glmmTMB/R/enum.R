## Auto generated - do not edit by hand
.valid_link <- c(
  log                 = 0,
  logit               = 1,
  probit              = 2,
  inverse             = 3,
  cloglog             = 4,
  identity            = 5,
  sqrt                = 6
)
.valid_family <- c(
  gaussian = 0,
  binomial = 100,
  betabinomial =101,
  beta =200,
  ordbeta = 201,
  Gamma =300,
  poisson =400,
  truncated_poisson =401,
  genpois =402,
  compois =403,
  truncated_genpois =404,
  truncated_compois =405,
  nbinom1 =500,
  nbinom2 =501,
  truncated_nbinom1 =502,
  truncated_nbinom2 =503,
  t =600,
  tweedie = 700,
  lognormal = 800
)
.valid_covstruct <- c(
  diag = 0,
  us   = 1,
  cs   = 2,
  ar1  = 3,
  ou   = 4,
  exp = 5,
  gau = 6,
  mat = 7,
  toep = 8,
  rr = 9,
  homdiag = 10
)
.valid_zipredictcode <- c(
  corrected = 0,
  uncorrected = 1,
  prob = 2,
  disp = 3
)
.valid_prior <- c(
  normal = 0,
  t = 1,
  cauchy = 2,
  gamma = 10,
  beta = 20,
  lkj = 30
)
.valid_vprior <- c(
  beta = 0,
  betazi = 1,
  betad = 2,
  theta = 10,
  thetazi = 20,
  psi = 30
)
