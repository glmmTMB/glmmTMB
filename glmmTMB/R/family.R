## FIXME: I would like to use the following function instead of repeating
## the pattern, but I'm worried that lazy evaluation of arguments will
## cause all kinds of trouble
family_factory <- function(default_link,family,variance) {
    f <- function(link=default_link) {
        r <- list(family=family,link=link,variance=variance)
        r <- c(r,make.link(link))
        class(r) <- "family"
        return(r)
    }
    return(f)
}

## suppress code warnings for nbinom2; can't use .Theta <- NULL trick here ...
globalVariables(".Theta") 

## attempt to guess whether calling function has been called from glm.fit ...
in_glm_fit <- function() {
    vars <- ls(envir=parent.frame(2))
    all(c("coefold","control","EMPTY","good","nvars") %in% vars)
}

make_family <- function(x,link) {
    x <- c(x,list(link=link),make.link(link))
    ## stubs for Effect.default/glm.fit
    if (is.null(x$aic)) {
        x <- c(x,list(aic=function(...) NA_real_))
    }
    if (is.null(x$initialize)) {
        ## should handle log-links adequately
        x <- c(x,list(initialize=expression({mustart <- y+0.1})))
    }
    if (is.null(x$dev.resids)) {
        ## can't return NA, glm.fit is unhappy
        x <- c(x,list(dev.resids=function(y,mu,wt)  {
                     rep(0,length(y))
                 }))
    }

    class(x) <- "family"
    return(x)
}

## even better (?) would be to have a standalone list including
## name, default link, variance function, (optionally) initialize
## for each family

##' Family functions for glmmTMB
##'
##'
##' 
##' @aliases family_glmmTMB
##' @param link (character) link function for the conditional mean ("log", "logit", "probit", "inverse", "cloglog", or "identity")
##' @return returns a list with (at least) components
##' \item{family}{length-1 character vector giving the family name}
##' \item{link}{length-1 character vector specifying the link function}
##' \item{variance}{a function of either 1 (mean) or 2 (mean and dispersion
##' parameter) arguments giving a value proportional to the
##' predicted variance (scaled by \code{sigma(.)})
##' }
##' @details
##' If specified, the dispersion model uses a log link. Denoting the dispersion parameter
##' as phi=exp(eta) (where eta is the linear predictor from the dispersion model)
##' and the predicted mean as mu:
##'  \describe{
##'      \item{gaussian}{(from base R): constant variance=phi}
##'      \item{Gamma}{(from base R) phi is the shape parameter, i.e variance=mu*phi}
##'      \item{nbinom2}{variance increases quadratically with the mean (Hardin & Hilbe 2007),
##' i.e. variance=mu*(1+mu/phi)}
##'      \item{nbinom1}{variance increases linearly with the mean (Hardin & Hilbe 2007),
##' i.e. variance=mu*(1+phi)}
##'      \item{compois}{is the Conway-Maxwell Poisson parameterized with the exact mean
##'           which differs from the COMPoissonReg package (Sellers & Lotze 2015)}
##'      \item{genpois}{is the generalized Poisson distribution}
##'      \item{beta}{follows the parameterization of Ferrari and Cribari-Neto (2004) and the \code{betareg} package,
##'     i.e. variance=mu*(1-mu)}
##' }
##' @references
##' \itemize{
##' \item Ferrari SLP, Cribari-Neto F (2004). "Beta Regression for Modelling Rates and Proportions." \emph{J. Appl. Stat.}  31(7), 799-815.
##' \item Hardin JW & Hilbe JM (2007). "Generalized linear models and extensions." Stata Press.
##' \item Sellers K & Lotze T (2015). "COMPoissonReg: Conway-Maxwell Poisson (COM-Poisson) Regression". R package version 0.3.5. https://CRAN.R-project.org/package=COMPoissonReg
##' }
##' @export
##' @importFrom stats make.link
nbinom2 <- function(link="log") {
    theta_errstr <- "theta (nbinom parameter) neither passed as an argument nor stored in enviroment"
    missing_theta <- "one" ## or "stop" or "na"
    r <- list(family="nbinom2",
              variance=function(mu, theta) {
                if (missing(theta)) {
                    ## look in environment
                    if (!exists(".Theta")) {
                        theta <- switch(missing_theta,
                                        one=1,
                                        na=NA_real_,
                                        stop=stop(theta_errstr))
                    } else {
                        theta <- .Theta
                    }
                }
                return(mu*(1+mu/theta))
               },  ## variance function
              ## full versions needed for effects::mer.to.glm
              ## (so we can evaluate a glm)
              initialize = expression({
                  if (any(y < 0)) 
                      stop("negative values not allowed for the negative binomial family")
                  n <- rep(1, nobs)
                  mustart <- y + (y == 0)/6
              }),
              dev.resids = function (y, mu, wt, theta)  {
        if (missing(theta)) {
            if (!exists(".Theta")) {
                theta <- switch(missing_theta,
                                na=NA_real_,
                                one=1,
                                stop=stop(theta_errstr))
            } else {
                theta <- .Theta
            }
        }
        return(2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta))))
    })
    return(make_family(r,link))
}

#' @rdname nbinom2
#' @export
nbinom1 <- function(link="log") {
    r <- list(family="nbinom1",
              variance=function(mu,alpha) {
        ## Effect stub (can't return 0 or NA or glm.fit will complain)
        ## FIXME: retrieve dispersion in environment?
        if (missing(alpha)) return(rep(1e-16,length(mu)))
        mu*(1+alpha)
    })
    return(make_family(r,link))
}

#' @rdname nbinom2
#' @export
compois <- function(link="log") {
    r <- list(family="compois",
           variance=function(mu,phi) {
               if (length(phi)==1) phi <- rep(phi, length=length(mu))
               .Call("compois_calc_var", mu, 1/phi, PACKAGE="glmmTMB")
          })
    return(make_family(r,link))
}

#' @rdname nbinom2
#' @export
truncated_compois <- function(link="log") {
    r <- list(family="truncated_compois",
           variance=function(mu,phi) {
             stop("variance for truncated compois family not yet implemented")
           })
    return(make_family(r,link))
}

#' @rdname nbinom2
#' @export
genpois <- function(link="log") {
    r <- list(family="genpois",
           variance=function(mu,phi) {
               mu*phi
           })
    return(make_family(r,link))
}

#' @rdname nbinom2
#' @export
truncated_genpois <- function(link="log") {
    r <- list(family="truncated_genpois",
           variance=function(mu,phi) {
             stop("variance for truncated genpois family not yet implemented")
          })
    return(make_family(r,link))
}

#' @rdname nbinom2
#' @export
truncated_poisson <- function(link="log") {
	r <- list(family="truncated_poisson",
           variance=function(lambda) {
           (lambda+lambda^2)/(1-exp(-lambda)) - lambda^2/((1-exp(-lambda))^2)
           })
        return(make_family(r,link))
}

#' @rdname nbinom2
#' @export	       	
truncated_nbinom2 <- function(link="log") {
    r <- list(family="truncated_nbinom2",
           variance=function(mu,theta) {
               stop("variance for truncated nbinom2 family not yet implemented")
         })
    return(make_family(r,link))
}

#' @rdname nbinom2
#' @export
truncated_nbinom1 <- function(link="log") {
    r <- list(family="truncated_nbinom1",
           variance=function(mu,alpha) {
               stop("variance for truncated nbinom1 family not yet implemented")
           })
    return(make_family(r,link))
}

## similar to mgcv::betar(), but simplified.
## variance has only one parameter; full variance is mu*(1-mu)/(1+phi) =
## sigma(.)*family(.)$variance(mu)
## initialize() tests for legal response values and sets trivial mustart
#' @rdname nbinom2
#' @export
beta_family <- function(link="logit") {
    ## note *internal* name must still be "beta",
    ## unless/until it's changed in src/glmmTMB.cpp (and R/enum.R is rebuilt)
    r <- list(family="beta",
              variance=function(mu) { mu*(1-mu) },
              initialize=expression({
                  if (any(y <= 0 | y >= 1)) 
                      stop("y values must be 0 < y < 1")
                  mustart <- y
              }))
    return(make_family(r,link))
}

## fixme: better name?

#' @rdname nbinom2
#' @export
betabinomial <- function(link="logit") {
    r <- list(family="betabinomial",
              variance=function(mu,phi) {
        stop("variance for betabinomial family not yet implemented")
    },
    initialize = binomial()$initialize)
    return(make_family(r,link))
}

#' @rdname nbinom2
#' @export
tweedie <- function(link="log") {
    r <- list(family="tweedie",
           variance=function(mu,phi,p) {
               stop("variance for tweedie family not yet implemented")
         })
    return(make_family(r,link))
}

## t not yet implemented in 
## t_family <- function(link="identity") {
##     ## FIXME: right now t behaves just like gaussian(); variance()
##     ## returns a value *proportional* to the variance
##     r <- list(family="t",link=link,
##               variance=function(mu) {
##         rep.int(1,length(mu))
##     })
## }

#' List model options that glmmTMB knows about
#'
#' @note these are all the options that are \emph{defined} internally; they have not necessarily all been \emph{implemented} (FIXME!)
#' @param what (character) which type of model structure to report on
#' ("all","family","link","covstruct")
#' @param check (logical) do brute-force checking to test whether families are really implemented (only available for \code{what="family"})
#' @return if \code{check==FALSE}, returns a vector of the names (or a list of name vectors) of allowable entries; if \code{check==TRUE}, returns a logical vector of working families
#' 
#' @export
getCapabilities <- function(what="all",check=FALSE) {
    if (!check) {
        switch(what,
               all=lapply(list(family=.valid_family,link=.valid_link,
                        covstruct=.valid_covstruct),names),
               family=names(.valid_family),
               link=names(.valid_link),
               covstruct=names(.valid_covstruct),
               stop(sprintf("unknown option %s",what)))
    } else {
        ## run dummy models to see if we get a family-not-implemented error
        if (what!="family") stop("'check' option only available for families")
        families <- names(.valid_family)
        family_OK <- setNames(rep(TRUE,length(.valid_family)),families)
        y <- 1:3 ## dummy
        for (f in families) {
            tt1 <- utils::capture.output(tt0 <- suppressMessages(suppressWarnings(
                try(glmmTMB(y~1,
                            family=list(family=f,link="identity")),
                    silent=TRUE))))
            family_OK[f] <- !(inherits(tt0,"try-error") &&
                              grepl("Family not implemented!",tt0))

        }
        return(family_OK)
    }
}
