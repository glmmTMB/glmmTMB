## we want to try to compute a reduced-rank neg log likelihood in two different ways

library(glmmTMB)
library(reformulas)
library(Matrix)
library(RTMB)

to_factormat <- function(theta, d) {
    nr <- (length(theta) + choose(d,2)) / d
    m <- matrix(0, nrow = nr, ncol = d)
    m[row(m) >= col(m)] <- theta
    m
}

from_factormat <- function(x) {
    x[row(x) >= col(x)]
}

## Log determinant example from RTMB
logdet <- ADjoint(
    function(x) determinant(x, log=TRUE)$modulus,
    function(x, y, dy) t(solve(x)) * dy,
    name = "logdet")


set.seed(101)
ngrp <- 30
nlev <- 30
d <- 3
ntheta <- nlev*d - choose(d,2)
nu <- d*nlev
nb <- ngrp*nlev
dd <- expand.grid(f = factor(1:ngrp),
                  g = factor(1:nlev))
form <- y ~ 1 + rr(0 + f | g, d = d)
dd$y <- simulate_new(RHSForm(form, as.form = TRUE),
                     newdata = dd, 
                     family = gaussian,
                     newparams = list(beta = 0,
                                      theta = rnorm(ntheta),
                                      betadisp = -1))[[1]]

rt <- mkReTrms(findbars(form), fr = model.frame (~ f + g, data = dd), calc.lambda = FALSE)
Z <- t(rt$Zt)
m1 <- glmmTMB(y ~ 1 + rr(f | g, d = d), 
              data = dd, 
              family = gaussian,
              control = glmmTMBControl(optCtrl = list(iter.max=1000,
                                                      eval.max=1000)))

## 1. using spherical random effects
f_spher <- function(par) {
    getAll(par, dd)
    Lambda <- matrix(0, nrow = nlev, ncol = d)
    Lambda[row(Lambda) >= col(Lambda)] <- theta
    b <- Lambda %*% matrix(u, nrow = d)
    mu <- beta0 + Z %*% c(b)
    sd1 <- exp(logsd)
    -sum(dnorm(y, drop(mu), sd1, log = TRUE)) - sum(dnorm(u, log = TRUE))
}
theta0 <- rnorm(d*nlev - choose(d,2))
par0 <- list(beta0 = 0, logsd = 0, theta = rnorm(ntheta),
                 u = rep(0, nu))
f_spher(par0)
gc()
ff0 <- MakeADFun(f_spher, par0)
ff0$fn()

## why does this go up?
ff1 <- MakeADFun(f_spher, par0, random = "u", silent = TRUE)
ff1$fn()


## 2. using non-spherical random effects
f_nonspher <- function(par) {
    getAll(par, dd)
    Phi <- to_factormat(theta, d)
    Phi[row(Phi) >= col(Phi)] <- theta
    mu <- drop(beta0 + Z %*% b)
    sd1 <- exp(logsd)
    u <- drop(Phi %*% matrix(b, nrow = d))
    logdetphi <- logdet(crossprod(Phi))
    -sum(dnorm(y, mu, sd1, log = TRUE)) + (sum(u^2) - logdetphi + d*log(2*pi))/2
}

par1 <- list(beta0 = 0, logsd = 0, theta = rnorm(ntheta),
             b = rep(0, nlev*ngrp))
f_nonspher(par1)
ff2 <- MakeADFun(f_nonspher, par1)
ff2$fn()
ff3 <- MakeADFun(f_nonspher, par1, random = "b", silent = TRUE)
ff3$fn()
## lots of memory, > ff2$fn() [-207.6635 -> 349.1435]

### playing with inverses

d <- 2
n <- 20
Lambda <- matrix(1:(d*n), ncol = d)
Sigma <- tcrossprod(Lambda)
Sigma_inv <- MASS::ginv(Sigma)
Phi <- MASS::ginv(Lambda)
all.equal(crossprod(Phi), Sigma_inv)

S <- svd(Lambda)
Phi_2 <- S$v[, 1:d] %*% ((1/S$d[1:d]) * t(S$u[, 1:d]))
all.equal(Phi, Phi_2)
