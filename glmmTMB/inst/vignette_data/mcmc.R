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

fm1 <- glmmTMB(count ~ mined + (1|site),
    zi=~mined,
    family=poisson, data=Salamanders)

## @knitr setup

## FIXME: is there a better way for user to extract full coefs?
rawcoef <- with(fm1$obj$env,last.par[-random])
names(rawcoef) <- make.names(names(rawcoef), unique=TRUE)
## log-likelihood function 
## (run_MCMC wants *positive* log-lik)
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

library(tmbstan)
t2 <- system.time(m2 <- tmbstan(fm1$obj, seed = 101))

## @knitr stanhacks

## functions to reduce the size of stored Stan-type objects
## obsolete?
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
rstan::traceplot(m2, pars=c("beta","betazi","theta"))
dev.off()

## @knitr tmbstan_pairsplot
png("tmbstan_pairsplot.png")
suppressWarnings(
    pairs(m2, pars = c("beta", "betazi"), gap = 0)
)
dev.off()

## @knitr diagnostics
dp2 <- bayestestR::diagnostic_posterior(m2)

## @knitr sleepstudy_tmbstan

data("sleepstudy", package = "lme4")
fm2 <- glmmTMB(Reaction ~ Days + (Days | Subject), data = sleepstudy)
t3 <- system.time(m3 <- tmbstan(fm2$obj, seed = 101))

## @knitr sleepstudy_diag
dp3 <- bayestestR::diagnostic_posterior(m3)

## @knitr sleepstudy_traceplot

png("sleepstudy_traceplot.png")
rstan::traceplot(m3, pars=c("beta","betadisp","theta"))
dev.off()

png("sleepstudy_traceplot_theta3.png")
rstan::traceplot(m3, pars=c("theta[3]"))
dev.off()

## @knitr sleepstudy_tmbstan_bounds

sdrsum <- summary(fm2$sdr)
par_est <- sdrsum[,"Estimate"]
par_sd <- sdrsum[,"Std. Error"]
t4 <- system.time(m4 <- tmbstan(fm2$obj,
                                lower = par_est - 5*par_sd,
                                upper = par_est + 5*par_sd,
                                seed = 101))

## @knitr sleepstudy_bounds_diag
dp4 <- bayestestR::diagnostic_posterior(m4)

## @knitr sleepstudy_bounds_traceplot
png("sleepstudy_traceplot_bounds.png")
rstan::traceplot(m4, pars=c("beta","betadisp","theta"))
dev.off()

png("sleepstudy_traceplot_bounds_theta3.png")
rstan::traceplot(m4, pars=c("theta[3]"))
dev.off()

## @knitr trans_param
samples4 <- as.data.frame(extract(m4, pars=c("beta","betadisp","theta")))
colnames(samples4) <- c(names(fixef(fm2)$cond),
                  "log(sigma)",
                  c("log(sd_Intercept)", "log(sd_Days)", "cor"))
samples4$cor <- sapply(samples4$cor, get_cor)

## @knitr save_all

## use version=2 to allow compatibility pre-3.5.0
## DON'T save tmbstan models (m[0-9]); even with size-hacking, not small enough.
## since PNG file is saved, we don't really need it
save(list = c("m1", ls(pattern = "(t|dp|samples)[0-9]+")), file="mcmc.rda", version=2)
