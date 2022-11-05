remotes::install_github("glmmTMB/glmmTMB/glmmTMB@denovo_simulation")
library(glmmTMB)

## OR (fork)/clone/pull glmmTMB repo; git checkout -b denovo_simulation
## devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")

## construct covariates for model (also see tools in `faux` package?)
dd <- expand.grid(trt = factor(c("A", "B")),
                  id = factor(1:10),
                  time = 1:6)
set.seed(101)
## response var needs to be something legal (e.g. in the right domain/
##  won't throw errors during initialization), but not necessarily sensible
dd$y <- rnbinom(nrow(dd), mu = 5, size = 2)

## construct data object required to build TMB object
r1 <- glmmTMB(y~trt*time+(1+time|id),
              data = dd,
              family = nbinom2,
              doFit = FALSE)
## construct TMB object, but don't fit it
r2 <- fitTMB(r1, doOptim = FALSE)
## inspect parameter vector
r2$env$par

## parameter vector is a bit odd: not a list, but a vector with
##  *duplicated* names corresponding to elements of the parameter lsit
set_pars <- function(base, ...) {
    ## FIXME: check for name matches, length matches etc.
    L <- list(...)
    for (nm in names(L)) {
        base[names(base) == nm] <- L[[nm]]
    }
    return(base)
}
## check fixed-effect predictor vars
colnames(X <- r1$data.tmb$X)
beta <- c(1, 2, 0.1, 0.2)
summary(exp(X %*% beta))
dim(Z <- r1$data.tmb$Z)

## would be convenient to have
##  a function for mapping relevant elements of theta to corr
##  for unstructured corr matrix, n > 2, this is the 
## inverse function of get_cor() -- not obvious!
## see covstruct vignete, "Mapping" section
## http://glmmtmb.github.io/glmmTMB/articles/covstruct.html

## for *2x2 matrices only* we only need to map the *last* element
## of theta to a single correlation parameter
## (again see covstruct vignette)
rho_to_theta <- function(rho) rho/sqrt(1-rho^2)
## test
stopifnot(all.equal(get_cor(rho_to_theta(-0.2)), -0.2))

## cor matrix -> cov matrix (for variances 1, 2)
D <- sqrt(diag(1:2))
Sigma <- D %*% matrix(c(1,-0.2, -0.2, 1), nrow = 2) %*% D
cov2cor(Sigma) ## double-check


## do want to generate this to test out expected range of X beta + Z b
b <- MASS::mvrnorm(1, mu = rep(0,20),
                   Sigma = Matrix::.bdiag(replicate(10,
                                                    Sigma,
                                                    simplify = FALSE)))

## unnecessary -- testing handmade MVnorm generation
cc <- chol(Sigma)
crossprod(cc)
## test
csims <- t(t(cc) %*% matrix(rnorm(20000), ncol = 10000))
var(csims)

## typical range of mean of fixed + RE on count scale:
## 0.0001 - 3000
summary(musim0 <- exp(drop(as.matrix((X %*% beta + Z %*% b)))))

## see ?family_glmmTMB for definitions of 'betad' -- usually
## log(dispersion parameter)
## we'll use phi = 1 -> betad = log(1) below
## i.e. geometric distribution
ysim0 <- rnbinom(length(musim0),
                 mu = musim0,
                 size = 1)

## range: 0 - 6000
summary(ysim0)

pp <- set_pars(r2$env$last.par, beta = c(1,2, 0.1, 0.2), betad = log(1),
               ## first elements of theta are on log-sd scale
               theta = c(log(sqrt(c(1,2))), rho_to_theta(-0.2)))

sim <- r2$simulate(pp)

## range is way too big (max = 4e9)
summary(sim$yobs)

## FIXME/TO DO:
## try this on simpler examples, e.g.
## - Gaussian (n.b. exp(betad) is residual **variance** ,
##    not sd, for gaussian family)
## - no random effects

