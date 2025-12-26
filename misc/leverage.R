## This file contains:
##
## 1. A simple standalone RTMB example that calulates the entire hat
##    matrix.
##
## 2. A general leverage function in glmmTMB based on RTMB with an
##    option to only calculate the diagonal in an efficent manner.
##
## NOTE: Assumes glmmTMB to be compiled in a way compatible with RTMB !
##   *  -DTMBAD_INDEX_TYPE=uint64_t
##   *  NO OPENMP !

## NOTES: this works, more or less, but for some reason the memory requirements
##  explode when we try it with a reduced-rank [ ... + rr(something) ...] model
##
## more tests on memory usage etc. can be found in misc/: see misc/memusg_master.R, misc/glmmTMB_mem.R
## (haven't looked at these in a long time! They use the `memprof` package to evaluate scaling etc.
## of memory requirements ...


######################################################################
## A simple standalone RTMB example
######################################################################

do_slow <- FALSE

library(RTMB)
library(Matrix)

parameters <- list(
    mua=0,          ## Mean slope
    sda=1,          ## Std of slopes
    mub=0,          ## Mean intercept
    sdb=1,          ## Std of intercepts
    sdeps=1,        ## Residual Std
    a=rep(0, 50),   ## Random slope by chick
    b=rep(0, 50)    ## Random intercept by chick
)
f <- function(parms) {
    getAll(ChickWeight, parms, warn=FALSE)
    ## Optional (enables extra RTMB features)
    weight <- OBS(weight)
    ## Initialize joint negative log likelihood
    nll <- 0
    ## Random slopes
    nll <- nll - sum(dnorm(a, mean=mua, sd=sda, log=TRUE))
    ## Random intercepts
    nll <- nll - sum(dnorm(b, mean=mub, sd=sdb, log=TRUE))
    ## Data
    predWeight <- a[Chick] * Time + b[Chick]
    nll <- nll - sum(dnorm(weight, predWeight, sd=sdeps, log=TRUE))
    ## Get predicted weight uncertainties
    ADREPORT(predWeight)
    ## Return
    nll
}
obj <- MakeADFun(f, parameters, random=c("a", "b"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
## Leverage calc
TMB::config(tmbad.atomic_sparse_log_determinant=0)
p <- obj$env$parList()
p$weight <- ChickWeight$weight
F <- GetTape(MakeADFun(f, p)) ## Objective: (theta, u, obs) -> nll
## Indices are as follows:
i_theta <- 1:5
i_u <- 6:105
i_obs <- 106:683
## Apply LA to integrate random effects out. Then opimize wrt theta
That <- F$laplace(i_u)$newton(i_theta) ## That: (obs) -> theta_hat
## Likewise, find Uhat: (theta, obs) -> u_hat
Uhat <- F$newton(i_u)
## Yhat tape
Yhat <- GetTape(MakeADFun(f, p, ADreport=TRUE)) ## (theta, u, obs) -> muhat
## Glue the stuff: (obs) -> muhat(theta(obs), u(obs), obs)
Hat <- MakeTape(function(y) {
    that <- That(y)
    uhat <- Uhat(c(that, y))
    Yhat(c(that, uhat, y))
}, p$weight)
## Entire hat matrix
zx <- Hat$jacobian(Hat$par())
## TODO: Compare with numDeriv::jacobian

L_rtmb <- diag(zx)

######################################################################
## glmmTMB implementation
######################################################################
library(glmmTMB)
library(RTMB)
library(Matrix)

## @param diag Get diagonal only?
leverage <- function(fm, diag=TRUE) {
    has.random <- any(fm$obj$env$lrandom())
    obj <- fm$obj
    ## We mess with these... (cleanup on exit!)
    restore.on.exit <- c("ADreport",
                         "parameters",
                         "data")
    oldvars <- sapply(restore.on.exit, get, envir=obj$env, simplify=FALSE)
    restore.oldvars <- function(){
        for(var in names(oldvars)) assign(var, oldvars[[var]], envir=obj$env)
    }
    on.exit({restore.oldvars(); obj$retape()})
    ## #################################################################
    ## Get derivatives of prediction
    ##
    ##    mu_hat( b_hat( theta_hat(yobs), yobs) , theta_hat(yobs) )
    ##
    ## wrt. yobs.
    ## Note the three partial derivative 'paths' to consider.
    ## We can get this derivative by
    ##  1. yobs -> theta_hat(yobs)
    ##  2. theta -> mu_hat( b_hat( theta, yobs) , theta ) [fixed yobs]
    ##  3. yobs -> mu_hat( b_hat( theta, yobs) , theta ) [fixed theta]
    ## #################################################################
    ##parhat <- obj$env$last.par.best
    parhat <- fm$fit$parfull
    pl <- obj$env$parList(par=parhat)
    yobs <- obj$env$data$yobs
    obj$env$parameters <- pl
    theta <- parhat[!obj$env$lrandom()]  ## ALL top-level parameters
    b <- parhat[obj$env$lrandom()]
    if (!is.null(obj$env$spHess)) {
        Hbb <- obj$env$spHess(parhat, random=TRUE) ## Needed later for RE models
    }
    ## #################################################################
    ## 1. Get partial derivatives of theta_hat wrt to yobs
    ## Note: length(yobs) much greater that length(theta)
    ##       ==> Reverse mode AD is suitable !
    ## #################################################################
    ## Move 'yobs' from data -> parameters (preserve C++ template order!)
    obj$env$parameters <- c(list(yobs = yobs), obj$env$parameters)
    obj$env$data$yobs <- NULL
    obj$retape()
    ## New objective parameters: (yobs, b, theta)
    nobs <- length(obj$env$parameters$yobs)
    nb <- length(obj$env$random)
    ntheta <- length(obj$env$par) - nobs - nb
    TMB::config(tmbad.atomic_sparse_log_determinant=0, DLL="RTMB") ## TMB FIXME
    F <- GetTape(obj)
    r <- obj$env$random ## Which are random
    p <- tail(1:(nobs+ntheta), ntheta) ## Which are parameters *after* removing random
    ThetaHat <- F$laplace(r)$newton(p)
    J <- ThetaHat$jacobian(ThetaHat$par())
    ## Extra stuff we need in (3)
    F. <- F$jacfun() ## (yobs, [b], theta) -> (yobs, [b], theta)
    F. <- MakeTape(function(y) F.( c(y, parhat) ) [r] , yobs) ## yobs -> b
    Hby <- F.$jacfun(sparse=TRUE)(yobs)
    ## #################################################################
    ## 2. Get partial derivatives of mu_hat wrt to theta for *fixed* yobs
    ## Note: length(mu) much greater that length(theta)
    ##       ==> Forward mode AD is suitable !
    ## #################################################################
    obj$env$data$yobs <- yobs
    obj$env$parameters$yobs <- NULL
    obj$retape()
    r <- obj$env$random ## Which are now random
    F <- GetTape(obj)
    Bhat <- F$newton(r) ## theta -> bhat
    obj$env$data$doPredict <- as.double(1) ## Enable prediction of 'mu'
    obj$env$data$whichPredict <- as.double(1:nobs)
    obj$env$ADreport <- TRUE ## Modify return value from Cpp
    obj$retape()
    F <- GetTape(obj) ## (b, theta) -> mu
    ## This doesn't work:
    ## MuHat <- MakeTape(function(theta)F(c(Bhat(theta), theta)), theta)
    MuHat <- MakeTape(function(theta) {
        par <- advector(nb + ntheta) ## glmmTMB mixes order of parameters and random effects...
        r <- obj$env$lrandom()
        par[r] <- Bhat(theta)
        par[!r] <- theta
        F(par)
    } , theta)
    ## 'Adjoint trick'
    T2 <- MakeTape(function(weight) {
        WT <- MakeTape(function(th) sum(MuHat(th) * weight), theta)
        WT$jacfun()(advector(theta))
    }, rep(1,nobs))
    J2 <- T2$jacobian(yobs)
    if (diag) {
        term1 <- colSums(J*J2)
    } else {
        term1 <- t(J) %*% J2
    }
    ## #################################################################
    ## 3. Get partial derivatives of mu_hat wrt yobs for fixed theta
    ## Note: Tricky!
    ##       
    ## #################################################################
    term2 <- 0
    if (has.random) {
        F2 <- MakeTape(function(b) {
            par <- advector(nb + ntheta) ## glmmTMB mixes order of parameters and random effects...
            r <- obj$env$lrandom()
            par[r] <- b
            par[!r] <- theta
            F(par)
        }, b) ## (b) -> mu
        Hmb <- F2$jacfun(sparse=TRUE)(b) ## sparse deriv mu wrt b
        ## Implicit function theorem gives final partial deriv matrix:
        ##   - Hby %*% solve(Hbb) %*% Hbm
        ## of which we need the diagonal.
        ## Because mu and yobs link to the same random effects, all required b-cliques are part of Hbb !
        ## It follows that we can replace solve(Hbb) by its subset iH !
        if (diag) {
            if (length(b) == 0) {
                iH <- solve(Hbb)
            } else {
                iH <- TMB:::solveSubset(Hbb)
            }
            browser()
            term2 <- -colSums( Hby * (  iH %*% t(Hmb) ) )
        } else {
            term2 <- -t(Hby) %*% solve(Hbb, t(Hmb))
        }
    }
    term1 + term2
}

## Test it
fm <- glmmTMB(weight ~ diag(Time | Chick) + Time, data=ChickWeight)
leverage(fm)
## Warning message:
## In GetTape(obj) : Permanently changing the global pointer of DLL 'glmmTMB'

## Check the first element:
## i <- 1
## myf <- function(eps) {
##     ChickWeight$weight[i] <- ChickWeight$weight[i] + eps
##     fm <- glmmTMB(weight ~ diag(Time | Chick) + Time, data=ChickWeight)
##     fitted(fm)
## }
## myf(0)
## numDeriv::jacobian(myf, 0)


fm0 <- glmmTMB(weight ~ Time, data=ChickWeight)
leverage(fm0)
gm0 <- lm(weight ~ Time, data = ChickWeight)
all.equal(leverage(fm0), unname(hatvalues(gm0)))

fm1 <- glmmTMB(weight ~ Time + diag(1 + Time|Chick), data=ChickWeight)
head(predict(fm1))
L_gt <- leverage(fm1)
## internals of model **NOT** restored ...
head(predict(fm0))
fm1 <- lme4::lmer(weight ~ Time + (1 + Time || Chick), data=ChickWeight)
L_l4 <- hatvalues(fm1)

plot(L_gt,L_l4)

if (do_slow) {
    ## ridiculous brute force (finite diffs)
    ## refit
    fm0 <- glmmTMB(weight ~ Time + diag(1 + Time|Chick), data=ChickWeight)
    nn <- nrow(ChickWeight)
    L_brute_gt <- numeric(nn)
    L_brute_l4 <- numeric(nn)
    eps <- 0.001
    if (interactive()) pb <- txtProgressBar(max=nn, style = 3)
    for (i in 1:nn) {
        if (interactive()) setTxtProgressBar(pb, i)
        newresp <- ChickWeight$weight
        newresp[i] <- newresp[i]+eps
        ## test fast_refit branch???
        L_brute_gt[i] <- (predict(refit(fm0, newresp))[i] - predict(fm0)[i])/eps
        L_brute_l4[i] <- (predict(refit(fm1, newresp))[i] - predict(fm1)[i])/eps
    }
    if (interactive()) close(pb)
    levmat <- cbind(L_brute_gt, L_brute_l4, L_gt, L_l4, L_rtmb)
    saveRDS(levmat, file = "levmat.rds")
} else {
    levmat <- readRDS("levmat.rds")
}
pairs(levmat, gap = 0,
      lower.panel = function(...) { points(...); abline(a=0, b=1, col = 2) },
      upper.panel = NULL)
