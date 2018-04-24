##' Compute likelihood profile confidence intervals of a TMB object by root-finding
##' (generalized from TMB::tmbprofile)
##' 
##' @param obj a fitted glmmTMB object
##' @param name parameter index/name
##' @param target desired deviation from minimum log-likelihood. Default
##' is set to retrieve the 95% likelihood profile confidence interval,
##' if the objective function is a negative log-likelihood function
##' @param lincomb linear combination of parameters
##' @param parm.range lower and upper limits; if \code{NA},
##' a value will be guessed based on the parameter value and \code{sd.range}
##' @param sd.range in the absence of explicit \code{parm.range} values,
##' the range chosen will be the parameter value plus or minus \code{sd.range}.
##' May be specified as a two-element vector for different ranges below and
##' above the parameter value.
##' @param trace report information?
##' @param continuation use continuation method, i.e. set starting parameters for non-focal parameters to solutions from previous fits?
##' @return a two-element numeric vector containing the lower and upper limits (or \code{NA} if the target is not achieved in the range), with an attribute giving the total number of function iterations used
##' @importFrom stats uniroot qchisq nlminb
##' @importFrom Matrix Diagonal
##' @importFrom TMB sdreport
##' @export
tmbroot <-
function (obj, name, target=0.5*qchisq(0.95,df=1),
          lincomb, parm.range = c(NA,NA),
          sd.range = 7,
          trace = FALSE,
          continuation = FALSE)
{
    ## continuation method works well for profiling, where
    ##  each fit starts "close" to previous values, but may be
    ##  counterproductive for root-finding, when we are jumping back
    ##  and forth ...
    restore.on.exit <- c("last.par.best", "random.start", "value.best", 
        "last.par", "inner.control", "tracemgc")
    oldvars <- sapply(restore.on.exit, get, envir = obj$env, 
        simplify = FALSE)
    restore.oldvars <- function() {
        for (var in names(oldvars)) assign(var, oldvars[[var]], 
            envir = obj$env)
    }
    on.exit(restore.oldvars())
    par <- obj$env$last.par.best
    if (!is.null(obj$env$random)) 
        par <- par[-obj$env$random]
    if (missing(lincomb)) {
        if (missing(name)) 
            stop("No 'name' or 'lincomb' specified")
        stopifnot(length(name) == 1)
        if (is.numeric(name)) {
            lincomb <- as.numeric(1:length(par) == name)
            name <- names(par)[name]
        }
        else if (is.character(name)) {
            if (sum(names(par) == name) != 1) 
                stop("'name' is not unique")
            lincomb <- as.numeric(names(par) == name)
        }
        else stop("Invalid name argument")
    }
    else {
        if (missing(name)) 
            name <- "parameter"
    }
    stopifnot(length(lincomb) == length(par))
    X <- Diagonal(length(lincomb))
    i <- which(lincomb != 0)[1]
    X[i, ] <- lincomb
    invX <- solve(X)
    direction <- invX[, i]
    C <- invX[, -i, drop = FALSE]
    that <- sum(lincomb * par)
    f <- function(x) {
        par <- par + x * direction
        newfn <- function(par0) {
            par <- par + as.vector(C %*% par0)
            obj$fn(par)
        }
        newgr <- function(par0) {
            par <- par + as.vector(C %*% par0)
            as.vector(obj$gr(par) %*% C)
        }
        obj$env$value.best <- Inf
        obj$env$inner.control$trace <- FALSE
        obj$env$tracemgc <- FALSE
        control <- list(step.min = 0.001)
        ans <- nlminb(start, newfn, newgr, control = control)
        if (continuation) start <<- ans$par
        conv <<- ans$convergence
        if (trace) 
            cat("Profile value:", ans$objective, "\n")
        ans$objective
    }
    f.original <- f
    f <- function(x) {
        y <- try(f.original(x), silent = TRUE)
        if (is(y, "try-error")) 
            y <- NA
        y
    }
    g <- function(x) {
        return(f(x)-v.0-target)
    }
    if (any(is.na(parm.range))) {
        sds <- TMB::sdreport(obj)
        sd0 <- drop(sqrt(lincomb %*% sds$cov.fixed %*% matrix(lincomb)))
        if (length(sd.range)==1) sd.range <- rep(sd.range,2)
        parm.range[is.na(parm.range)] <-
            c(-1,1)*sd0*sd.range[is.na(parm.range)]
    }
    ## need to set start in order for f() to work ...
    ## FIXME: check convergence code ...
    conv <- 0
    start <- rep(0, length(par) - 1)
    v.0 <- f(0) ## need to set v.0 for g() ...
    lwr.x <- g(parm.range[1])
    if (is.na(lwr.x) || lwr.x<0) {
        lwr <- list(root=NA,iter=0)
    } else {
        lwr <- uniroot(g,interval=c(parm.range[1],0))
    }
    ## reset for upper root-finding
    restore.oldvars()
    start <- rep(0, length(par) - 1)
    upr.x <- g(parm.range[2])
    if (is.na(upr.x) || upr.x<0) {
        upr <- list(root=NA,iter=0)
    } else {
        upr <- uniroot(g,interval=c(0,parm.range[2]))
    }
    ans <- c(lwr=that+lwr$root,upr=that+upr$root)
    attr(ans,"iter") <- lwr$iter+upr$iter
    return(ans)
}
