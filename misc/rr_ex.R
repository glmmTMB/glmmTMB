## we want to try to compute a reduced-rank neg log likelihood in two
## different ways:

## * with spherical random effects, parameterized by factor loading
## matrix (eta = X*beta + Z*Lambda*u, u ~ dnorm())

## * with nonspherical random effects, param by inverse factor loading
## matrix (eta = X*beta + Z*b, b ~ dmvnorm(Phi)

## The hope is that using nonspherical random effects will speed up
## leverage calculations by making the latent-variable Hessian
## sparser, although in preliminary tests below it seems to slow down
## the fitting machinery itself.

## Current status: still fighting with how/whether to subtract
## log(det(tcrossprod(Phi))) in the non-spherical log-likelihood ...
## I think I need it as a Jacobian term to account for moving from b =
## (Lambda u) to u (i.e., left-multiplying by Phi) ?  Still a bit
## confused; leave it out if we just want to match NLLs between
## spherical and non-spherical parameterizations, but need to include
## it in the random-effects/Laplace-approximated version?  At the
## moment the non-spherical stuff doesn't seem to be working; I
## believe the precision is going to zero/RE variance to infinity so
## that the latent variables can match the data more or less exactly
## ...

## u = vector of spherical latent vars
## n = length(u)
## b

## penalty for spherical random effects:
## NLL(u ~ MVN(0, I)) = -sum(dnorm(u)) = sum(u^2)/2 + n/2*log(2*pi)
##
## penalty for non-spherical random effects:
##   NLL(b ~ MVN(0, Sigma)) = -sum(dmvnorm(b[,i], Sigma0)) =
## sum(b[,i] %*% invSigma %*% t(b[,i])/2) + n*log(det(invSigma))
##   = sum(tcrossprod(b[,i] %*% Phi)) + ngrp*log(det(Phi))/2 + n*log(2*pi)/2

## TO DO:
## * check Hessian patterns (how?)
## * try to figure out why Laplace approx is giving different answers for
##   spherical vs nonspherical; what can we do with TMB objects?

library(glmmTMB)
library(reformulas)
library(Matrix)
library(RTMB)
library(rbenchmark)
library(vegan)  ## for procrustes()
do_slow <- FALSE

#' convert parameter vector to constrained factor loading matrix
#' @param theta parameter vector
#' @param d dimension (rank)
#' @examples
#' to_factormat(1:7, d = 2)
to_factormat <- function(theta, d) {
    nd <- (length(theta) + choose(d,2)) / d
    m <- matrix(0, nrow = nd, ncol = d)
    m[row(m) >= col(m)] <- theta
    m
}

#' convert constrained factor matrix to parameter vector
from_factormat <- function(x) {
    x[row(x) >= col(x)]
}

## Log determinant example from RTMB
logdet <- ADjoint(
    function(x) determinant(x, log=TRUE)$modulus,
    function(x, y, dy) t(solve(x)) * dy,
    name = "logdet"
)

set.seed(101)
## WARNING: we'll assume ngrp/nlev/d are fixed throughout this
## example, modifying them will break stuff
ngrp <- 23
nlev <- 30
d <- 3
ntheta <- nlev*d - choose(d,2)
nu <- d*ngrp     ## number of spherical LVs
nb <- nlev*ngrp  ## number of nonspher LVs
dd <- expand.grid(f = factor(1:nlev),
                  g = factor(1:ngrp))
form <- y ~ 1 + rr(0 + f | g, d = d)
dd$y <- simulate_new(RHSForm(form, as.form = TRUE),
                     newdata = dd, 
                     family = gaussian,
                     newparams = list(beta = 0,
                                      theta = rnorm(ntheta),
                                      betadisp = -1))[[1]]

fitfun <- function(dd) {
    glmmTMB(form,
            data = dd,
            family = gaussian,
            control = glmmTMBControl(optCtrl = list(iter.max=1000,
                                                    eval.max=1000)))
}

m1 <- fitfun(dd)

## construct Z matrix
rt <- mkReTrms(findbars(form), fr = model.frame (~ f + g, data = dd), calc.lambda = FALSE)
Z <- t(rt$Zt)

## 1. using spherical random effects
f_spher <- function(par) {
    getAll(par, dd)
    Lambda <- to_factormat(theta, d)
    b <- Lambda %*% matrix(u, nrow = d)
    mu <- drop(beta0 + Z %*% c(b))
    REPORT(mu)
    REPORT(b)
    nll <- -sum(dnorm(y, mu, exp(logsd), log = TRUE))
    nllpen <- -sum(dnorm(u, log = TRUE))
    nll + nllpen
}

## containerized version of MakeADFun
mk_f_spher <- function(par, dd, d, random = "u") {
    MakeADFun(f_spher, par, random = random, silent = TRUE)
}

theta0 <- rnorm(ntheta)
u0 <- rnorm(nu)
par0 <- list(beta0 = 0, logsd = 0, theta = theta0, u = u0)
nll0 <- f_spher(par0)
ff0 <- MakeADFun(f_spher, par0)
stopifnot(all.equal(nll0, ff0$fn()))

ff1 <- mk_f_spher(par0, dd, d)
ff1$fn()

pp1 <- ff1$env$parList()
pp1_b <- c(ff1$report()$b)

## 2. using non-spherical random effects
f_nonspher <- function(par) {
    getAll(par, dd)
    ## construct factor matrix with or without constraints?
    if (!raw_phi) {
        Phi <- t(to_factormat(theta, d))   ## zero-triangle
    } else {
        Phi <- t(matrix(theta, ncol = d))  ## no
    }
    mu <- drop(beta0 + Z %*% b)
    b <- matrix(b, nrow = nlev)
    u <- drop(Phi %*% b)
    nll <- -sum(dnorm(y, mu, exp(logsd), log = TRUE))
    nllpen <- -sum(dnorm(u, log = TRUE))
    if (jac_corr) {
        logdetphi <- ngrp*logdet(tcrossprod(Phi))
        nllpen <- nllpen - logdetphi  ## sign ???
    }
    REPORT(mu)
    REPORT(b)
    nll + nllpen 
}

mk_f_nonspher <- function(par, dd, d, random = "b", raw_phi = FALSE, jac_corr = TRUE) {
    MakeADFun(f_nonspher, par, random = random, silent = TRUE)
}

## check transformations, getting ready for non-spherical equivalents

## Phi0 %*% Lambda0 = I
Lambda0 <- to_factormat(theta0, d)
Phi0 <- MASS::ginv(Lambda0)
stopifnot(all.equal(Phi0 %*% Lambda0, diag(3)))

## Phi %*% b == u
b0 <- Lambda0 %*% matrix(u0, nrow = d)
stopifnot(all.equal(Phi0 %*% b0, matrix(u0, nrow = d)))

## for this section, keep full inverse factor loading matrix (without zeroing the upper triangle),
## to ensure comparability between original and inverse parameterizations
raw_phi <- TRUE
jac_corr <- FALSE
par1 <- list(beta0 = 0, logsd = 0, theta = t(Phi0), b = c(b0))
nll1 <- f_nonspher(par1)
## says we *DON'T* need Jacobian correction (I can talk myself into this),
## as we're still evaluating the log-likelihood of u|theta (not b|theta)
stopifnot(all.equal(nll1, nll0))

ff2 <- mk_f_nonspher(par1, dd, d, random = character(0), raw_phi = TRUE, jac_corr = FALSE)
ff2$fn()
stopifnot(all.equal(nll1, ff2$fn()))
pp2 <- ff2$env$parList()

## with random effects
ff3 <- mk_f_nonspher(par1, dd, d, raw_phi = TRUE, jac_corr = TRUE)
ff3$fn()
pp3 <- ff3$env$parList()

## check results (estimated b/u values)
## shouldn't these be identical? very weakly correlated
head(pp3$b)
head(pp1_b)
plot(ff3$report()$mu, ff1$report()$mu)

## nonspher much slower than spher for a single function evaluation
benchmark(spher = ff1$fn(), nonspher = ff3$fn())

## not as much difference for gradients, but still slower (9x)
benchmark(spher = ff1$gr(), nonspher = ff3$gr())

## now use identifiability constraint
raw_phi <- FALSE 
par1B <- list(beta0 = 0, logsd = 0, theta = from_factormat(t(Phi0)), b = c(b0))
f_nonspher(par1B)
## not very different from example with inverse factor loading
## matrix complete (i.e., not zeroing out the triangle: 3027.73 vs 3027.744
ff4 <- mk_f_nonspher(par1B, dd, d, raw_phi = FALSE, jac_corr = FALSE)
ff4$fn()  ## 92.01678 (vs 1854 for original parameterization with inner optimization ...)

## do optimization
ctrl <- list(iter.max = 1000, eval.max = 1000)
system.time(
    fit1 <- with(ff1, nlminb(par, fn, gr, control = ctrl))
)
system.time(
    fit4 <- with(ff4, nlminb(par, fn, gr))
)
## takes longer, but fits data much better

## RTMB agrees with glmmTMB fit
c(fit1$object, c(-1*logLik(m1)))

flatten <- function(x) c(as.matrix(x))
pmat <- cbind(nonspher=flatten(ff1$report()$mu),
              spher=flatten(ff4$report()$mu),
              glmmTMB = predict(m1), true = dd$y)
pairs(pmat, gap = 0)
hist(pmat[,"spher"]-pmat[,"true"])
## precision-matrix version is getting to the right place ...

## compare estimated factor loadings

simfun <- function() {
    dd_sim <- dd
    dd_sim$y <- simulate_new(RHSForm(form, as.form = TRUE),
                             newdata = dd, 
                             family = gaussian,
                             newparams = list(beta = 0,
                                              theta = rnorm(ntheta),
                                              betadisp = -1))[[1]]
    return(dd_sim)
}

fn <- "rr_sim_results.rds"
if (do_slow) {
    set.seed(101)
    nsims <- 50
    res <- list()
    for (i in 1:nsims) {
        dd_sim <- simfun()
        cat(i, "\n")
        th0 <- rnorm(d*nlev - choose(d,2))
        u0 <- rnorm(nu)
        LL <- to_factormat(th0, d)
        Phi0 <- MASS::ginv(Lambda0)
        b0 <- LL %*% matrix(u0, nrow = d)
        rpar0 <- list(beta0 = 0, logsd = 0, theta = th0, u = u0)
        rpar1 <- list(beta0 = 0, logsd = 0, theta = from_factormat(t(Phi0)), b = c(b0))
        trpar <- function(x) unlist(x[!names(x) %in% c("u", "b")])
        ff1sim <- mk_f_spher(rpar0, dd_sim, d)
        ff4sim <- mk_f_nonspher(rpar1, dd_sim, d)
        slik0 <- ff1sim$fn()
        slik1 <- ff4sim$fn()
        t0 <- system.time(rfit0 <- with(ff1, nlminb(par, fn, gr, control = ctrl)))
        t1 <- system.time(rfit1 <- with(ff4, nlminb(par, fn, gr, control = ctrl)))
        t2 <- system.time(rfit2 <- fitfun(dd_sim))
        res[[i]] <- tibble::lst(rpar0, rpar1, rfit0, rfit1, rfit2, slik0, slik1, t0, t1, t2)
        saveRDS(res, file = fn)
    }
} else {
    if (file.exists(fn)) res <- readRDS(fn)
}

## NEXT:
## * Compare with built-in rr()? (DONE)
## * Try setting starting values based on Dunn-Smyth glm resids?
## * Simulate new data each time?
## * see what happens with DEoptim? (i.e. is the likelihood surface horrible/multimodal?)

### playing with inverses

d0 <- 2  ## avoid modifying d when playing interactively 
n <- 20
Lambda <- matrix(1:(d*n), ncol = d0)
Sigma <- tcrossprod(Lambda)
Sigma_inv <- MASS::ginv(Sigma)
Phi <- MASS::ginv(Lambda)
all.equal(crossprod(Phi), Sigma_inv)

S <- svd(Lambda)
Phi_2 <- S$v[, 1:d0] %*% ((1/S$d[1:d0]) * t(S$u[, 1:d0]))
all.equal(Phi, Phi_2)

### does spider example work with new rr_inv covstruct (on branch??)
sppTot <- sort(
    tapply(spider_long$abund, list(spider_long$Species), sum),
    decreasing = TRUE
)
spiderDat_common <- subset(spider_long, Species %in% names(sppTot)[1:4])
spider_p1 <- glmmTMB(abund ~ Species + rr(Species + 0|id, d = 1),
                     family = poisson,
                     data=spiderDat_common)

spider_p2 <- glmmTMB(abund ~ Species + rr_inv(Species + 0|id, d = 1),
                     family = poisson,
                     data=spiderDat_common)

spider_p1$obj$fn()
spider_p2$obj$fn()
