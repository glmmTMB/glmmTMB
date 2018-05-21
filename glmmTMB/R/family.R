## FIXME: I would like to use the following function instead of repeating
## the pattern, but I'm worried that lazy evaluation of arguments will
## cause all kinds of trouble
family_factory <- function(default_link,family,variance) {
    f <- function(link=default_link) {
        r <- list(family=family,link=link,variance=variance)
        r <- c(r,make.link(link))
        return(r)
    }
    return(f)
}
## even better (?) would be to have a standalone list including
## name, default link, variance function, (optionally) initialize
## for each family

##' Family functions for glmmTMB
##' 
##' @aliases family_glmmTMB
##' @param link (character) link function for the conditional mean ("log", "logit", "probit", "inverse", "cloglog", or "identity")
##' @return returns a list with (at least) components
##' \item{family}{length-1 character vector giving the family name}
##' \item{link}{length-1 character vector specifying the link function}
##' \item{variance}{a function of either 1 (mean) or 2 (mean and dispersion
##' parameter) arguments giving the predicted variance}
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
##'     i.e. variance=mu*(1-mu)/(1+phi)}
##' }
##' @references
##' \itemize{
##' \item Ferrari SLP, Cribari-Neto F (2004). "Beta Regression for Modelling Rates and Proportions." \emph{J. Appl. Stat.}  31(7), 799-815.
##' \item Hardin JW & Hilbe JM (2007). "Generalized linear models and extensions." Stata Press.
##' \item Sellers K & Lotze T (2015). "COMPoissonReg: Conway-Maxwell Poisson (COM-Poisson) Regression". R package version 0.3.5. https://CRAN.R-project.org/package=COMPoissonReg
##' }
##' @importFrom stats make.link
##' @export
nbinom2 <- function(link="log") {
   r <- list(family="nbinom2",link=link,
           variance=function(mu,theta) {
               mu*(1+mu/theta)
       })
       return(c(r,make.link(link)))
}

#' @rdname nbinom2
#' @export
nbinom1 <- function(link="log") {
    r <- list(family="nbinom1",link=link,
              variance=function(mu,alpha) {
                  mu*(1+alpha)
              })
    return(c(r,make.link(link)))
}

#' @rdname nbinom2
#' @export
compois <- function(link="log") {
    r <- list(family="compois",link=link,
           variance=function(mu,phi) {
               if (length(phi)==1) phi <- rep(phi, length=length(mu))
               .Call("compois_calc_var", mu, 1/phi, PACKAGE="glmmTMB")
          })
    return(c(r,make.link(link)))
}

#' @rdname nbinom2
#' @export
truncated_compois <- function(link="log") {
    r <- list(family="truncated_compois",link=link,
           variance=function(mu,phi) {
             stop("variance for truncated compois family not yet implemented")
           })
    return(c(r,make.link(link)))
}

#' @rdname nbinom2
#' @export
genpois <- function(link="log") {
    r <- list(family="genpois",link=link,
           variance=function(mu,phi) {
               mu*phi
           })
    return(c(r,make.link(link)))
}

#' @rdname nbinom2
#' @export
truncated_genpois <- function(link="log") {
    r <- list(family="truncated_genpois",link=link,
           variance=function(mu,phi) {
             stop("variance for truncated genpois family not yet implemented")
          })
    return(c(r,make.link(link)))
}

#' @rdname nbinom2
#' @export
truncated_poisson <- function(link="log") {
	r <- list(family="truncated_poisson", link=link,
           variance=function(lambda) {
           (lambda+lambda^2)/(1-exp(-lambda)) - lambda^2/((1-exp(-lambda))^2)
           })
        return(c(r,make.link(link)))
}

#' @rdname nbinom2
#' @export	       	
truncated_nbinom2 <- function(link="log") {
    r <- list(family="truncated_nbinom2",link=link,
           variance=function(mu,theta) {
               stop("variance for truncated nbinom2 family not yet implemented")
         })
    return(c(r,make.link(link)))
}

#' @rdname nbinom2
#' @export
truncated_nbinom1 <- function(link="log") {
    r <- list(family="truncated_nbinom1",link=link,
           variance=function(mu,alpha) {
               stop("variance for truncated nbinom1 family not yet implemented")
           })
    return(c(r,make.link(link)))
}

## similar to mgcv::betar(), but simplified (variance has two parameters
##  rather than retrieving a variable from the environment); initialize()
##  tests for legal response values
#' @rdname nbinom2
#' @export
beta_family <- function(link="logit") {
    ## note *internal* name must still be "beta",
    ## unless/until it's changed in src/glmmTMB.cpp (and R/enum.R is rebuilt)
    r <- list(family="beta",link=link,
                variance=function(mu,phi) {
                    mu*(1-mu)/(1+phi)
                },
                initialize=expression({
                    if (any(y <= 0 | y >= 1)) 
                        stop("y values must be 0 < y < 1")
                }))
    return(c(r,make.link(link)))
}
## fixme: better name?

#' @rdname nbinom2
#' @export
betabinomial <- function(link="logit") {
    r <- list(family="betabinomial",
                link=link,
                variance=function(mu,phi) {
        stop("variance for betabinomial family not yet implemented")
    },
    initialize = binomial()$initialize)
    return(c(r,make.link(link)))
}

#' @rdname nbinom2
#' @export
tweedie <- function(link="log") {
    r <- list(family="tweedie",link=link,
           variance=function(mu,phi,p) {
               stop("variance for tweedie family not yet implemented")
         })
    return(c(r,make.link(link)))
}

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
