## @knitr run_MCMC

##' @param start starting value
##' @param V variance-covariance matrix of MVN candidate distribution
##' @param iterations total iterations
##' @param nsamp number of samples to store
##' @param burnin number of initial samples to discard
##' @param thin thinning interval
##' @param tune tuning parameters; expand/contract V
##' @param seed random-number seed
run_MCMC <- function(start,
                     V,   
                     logpost_fun,
                     iterations = 10000,
                     nsamp = 1000,
                     burnin = iterations/2,
                     thin = floor((iterations-burnin)/nsamp),
                     tune = NULL,
                     seed=NULL
                     ) {
    ## initialize
    if (!is.null(seed)) set.seed(seed)
    if (!is.null(tune)) {
        tunesq <- if (length(tune)==1) tune^2 else outer(tune,tune)
        V <-  V*tunesq
    }
    chain <- matrix(NA, nsamp+1, length(start))
    chain[1,] <- cur_par <- start
    postval <- logpost_fun(cur_par)
    j <- 1
    for (i in 1:iterations) {
        proposal = MASS::mvrnorm(1,mu=cur_par, Sigma=V)
        newpostval <- logpost_fun(proposal)
        accept_prob <- exp(newpostval - postval)
        if (runif(1) < accept_prob) {
            cur_par <- proposal
            postval <- newpostval
        }
        if ((i>burnin) && (i %% thin == 1)) {
            chain[j+1,] <- cur_par
            j <- j + 1
        }
    }
    chain <- na.omit(chain)
    colnames(chain) <- names(start)
    chain <- coda::mcmc(chain)
    return(chain)
}

## @knitr pkgs

library(glmmTMB)
library(tmbstan)
library(coda)

## @knitr fit1

data("sleepstudy",package="lme4")
fm1 <- glmmTMB(Reaction ~ Days + (Days|Subject),
               sleepstudy)


## @knitr setup

## FIXME: is there a better way for user to extract full coefs?
rawcoef <- with(fm1$obj$env,last.par[-random])
names(rawcoef) <- make.names(names(rawcoef),unique=TRUE)
## log-likelihood function 
## (MCMCmetrop1R wants *positive* log-lik)
logpost_fun <- function(x) -fm1$obj$fn(x)
## check definitions
stopifnot(all.equal(c(logpost_fun(rawcoef)),
                    c(logLik(fm1)),
          tolerance=1e-7))
V <- vcov(fm1,full=TRUE)

## @knitr do_run_MCMC
t1 <- system.time(m1 <- run_MCMC(start=rawcoef,
                                 V=V, logpost_fun=logpost_fun,
                                 seed=1001))

## @knitr do_tmbstan

## install.packages("tmbstan")
library(tmbstan)
t2 <- system.time(m2 <- tmbstan(fm1$obj))

## @knitr stanhacks

## functions to reduce the size of stored Stan-type objects
hack_size <- function(x, ...) {
    UseMethod("hack_size")
}

hack_size.stanfit <- function(x) {
    x@stanmodel <- structure(numeric(0), class="stanmodel")
    x@.MISC <- new.env()
    return(x)
}

hack_size.brmsfit <- function(x) {
    x$fit <- hack_size(x$fit)
    return(x)
}

hack_size.stanreg <- function(x) {
    x$stanfit <- hack_size(x$stanfit)
    return(x)
}
m2 <- hack_size(m2)

## @knitr tmbstan_traceplot

png("tmbstan_traceplot.png")
rstan::traceplot(m2, pars=c("beta","betad","theta"))
dev.off()

## @knitr save_all

## use version=2 to allow compatibility pre-3.5.0
## DON'T save m2; even with size-hacking, not small enough.
## since PNG file is saved, we don't really need it
save("m1","t1","t2", file="mcmc.rda", version=2)
