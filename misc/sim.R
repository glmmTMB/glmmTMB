library(glmmTMB)
library(lme4)

## OR (fork)/clone/pull glmmTMB repo; git checkout -b denovo_simulation
## devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")

##  a function for mapping relevant elements of theta to corr
##  for unstructured corr matrix, n > 2, this is the 
## inverse function of get_cor() 
## see covstruct vignete, "Mapping" section
## http://glmmtmb.github.io/glmmTMB/articles/covstruct.html

## C is a correlation matrix
## convert correlation matrix to scaled correlation parameters
put_cor <- function(C) {
    cc <- chol(C)
    cc2 <- t(cc %*% diag(1/diag(cc)))
    cc2[lower.tri(cc2)]
}

## testing put_cor
C <- matrix(c(1,  0.2,  0.1,
              0.2,  1, -0.2,
              0.1,-0.2,   1),
            3, 3)

## test: round-trip (almostl results in lower triangle only)
stopifnot(all.equal(get_cor(put_cor(C)),
                    C[lower.tri(C)]))

## for *2x2 matrices only* we only need to map the *last* element
## of theta to a single correlation parameter
## (again see covstruct vignette)
rho_to_theta <- function(rho) rho/sqrt(1-rho^2)
## test
stopifnot(all.equal(get_cor(rho_to_theta(-0.2)), -0.2))

## check rho_to_theta (specific to 2x2) vs more general solution

stopifnot(all.equal(rho_to_theta(-0.2), put_cor(matrix(c(1,-0.2,-0.2,1), 2))))

## parameter vector is a bit odd: not a list, but a vector with
##  *duplicated* names corresponding to elements of the parameter lsit
make_pars <- function(pars, ...) {
    ## FIXME: check for name matches, length matches etc.
    L <- list(...)
    for (nm in names(L)) {
        pars[names(pars) == nm] <- L[[nm]]
    }
    return(pars)
}

simfun <- function(form, data, pars, show_pars = FALSE, ...) {
    ## for now assume response variable is in data
    r1 <- glmmTMB(form,
              data = data,
              ...,
              doFit = FALSE)
## construct TMB object, but don't fit it
    r2 <- fitTMB(r1, doOptim = FALSE)
    if (show_pars) return(r2$env$last.par)
    pars <- do.call("make_pars",
                    c(list(r2$env$last.par), pars))
    r2$simulate(par = pars)$yobs
}


## example 1

data("sleepstudy", package = "lme4")
set.seed(101)
ss <-
    data.frame(sleepstudy,
               s1 = simfun(Reaction ~ Days, sleepstudy, 
                           pars = list(beta = c(250, 10), betad = 2*log(47))))

library(ggplot2)
ggplot(ss, aes(x = Days, colour = Subject)) +
    geom_line(aes(y=Reaction)) +
    geom_line(aes(y=s1), lty = 2)

## example 2

## set up metadata (one 'taxon' only)

## construct covariates for model (also see tools in `faux` package?)
dd <- expand.grid(trt = factor(c("A", "B")),
                  id = factor(1:10),
                  time = 1:6)
## response var needs to exist and be legal (i.e. within the domain
##  of the response distribution, probably not all zero)
##  but not necessarily sensible
dd$y <- rep(1, nrow(dd))

## parameter values (in human-understandable terms)
## make coef var approximately = 1 across groups (SDs the same as
##   the average mean and slope across control/treatment)
sdvec <-  sqrt(c(1.5,0.15))
corval <- -0.2

## first elements of theta
##  (corresponding to sds) are on log-sd scale; last elements describe
##  the correlation
thetavec <- c(log(sdvec), rho_to_theta(corval))

par1 <- list(beta = c(1,2, 0.1, 0.2),
             betad = log(1),  ## log(theta)
             theta = thetavec)
form1 <- y~trt*time+(1+time|id)


set.seed(101)
dd$y <- simfun(form1,
               data = dd,
               family = nbinom2,
               pars = par1)

range(dd$y)

## now reconstruct by hand (with some help from lme4 machinery)
set.seed(101)
## corresponding theta values 

X <- model.matrix(~trt*time, data  = dd)

## generate random effects values
rt <- mkReTrms(findbars(form1),
               model.frame(subbars(form1), data = dd))
Z <- t(rt$Zt)
C2 <- diag(sdvec) %*% matrix(c(1, corval, corval, 1), 2) %*% diag(sdvec)
lmer_thetavec <- t(chol(C2))[c(1,2,4)]

rt$Lambdat@x <- lmer_thetavec[rt$Lind]
u <- rnorm(nrow(rt$Lambdat))
b <- t(rt$Lambdat) %*% u

eta <- drop(X %*% par1$beta + t(rt$Zt) %*% b)
mu <- exp(eta)
y <- rnbinom(nrow(dd), size = 1, mu = mu)
range(y) ## range varies a lot


## JUNK/alternatives

## do want to generate this to test out expected range of X beta + Z b
b <- MASS::mvrnorm(1, mu = rep(0,20),
                   Sigma = Matrix::.bdiag(replicate(10,
                                                    Sigma,
                                                    simplify = FALSE)))


## FIXME/TO DO:
## try this on simpler examples, e.g.
## - Gaussian (n.b. exp(betad) is residual **variance** ,
##    not sd, for gaussian family)
## - Poisson?
## - no random effects
## - check against cases that lme4 simulate() can do
##   (may not get identical answers but should be identical
##   on average, reasonable ranges)


