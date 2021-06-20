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
utils::globalVariables(".Theta") 

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
##' @param link (character) link function for the conditional mean ("log", "logit", "probit", "inverse", "cloglog", "identity", or "sqrt")
##' @return returns a list with (at least) components
##' \item{family}{length-1 character vector giving the family name}
##' \item{link}{length-1 character vector specifying the link function}
##' \item{variance}{a function of either 1 (mean) or 2 (mean and dispersion
##' parameter) arguments giving a value proportional to the
##' predicted variance (scaled by \code{sigma(.)})
##' }
##' @details
##' If specified, the dispersion model uses a log link. Denoting the variance as \eqn{V}, the dispersion parameter
##' as \eqn{\phi=\exp(\eta)}{phi=exp(eta)} (where \eqn{\eta}{eta} is the linear predictor from the dispersion model), 
##' and the predicted mean as \eqn{\mu}{mu}:
##'  \describe{
##'      \item{gaussian}{(from base R): constant \eqn{V=\phi}{V=phi}}
##'      \item{Gamma}{(from base R) phi is the shape parameter. \eqn{V=\mu\phi}{V=mu*phi}}
##'       \item{ziGamma}{a modified version of \code{Gamma} that skips checks for zero values, allowing it to be used to fit hurdle-Gamma models}
##'      \item{nbinom2}{Negative binomial distribution: quadratic parameterization (Hardin & Hilbe 2007). \eqn{V=\mu(1+\mu/\phi) = \mu+\mu^2/\phi}{V=mu*(1+mu/phi) = mu+mu^2/phi}.}
##'      \item{nbinom1}{Negative binomial distribution: linear parameterization (Hardin & Hilbe 2007). \eqn{V=\mu(1+\phi)}{V=mu*(1+phi)}}
##'      \item{truncated_nbinom2}{Zero-truncated version of nbinom2: variance expression from Shonkwiler 2016. Simulation code (for this and the other truncated count distributions) is taken from C. Geyer's functions in the \code{aster} package; the algorithms are described in \href{https://cran.r-project.org/package=aster/vignettes/trunc.pdf}{this vignette}.}
##'      \item{compois}{Conway-Maxwell Poisson distribution: parameterized with the exact mean (Huang 2017), which differs from the parameterization used in the \pkg{COMPoissonReg} package (Sellers & Shmueli 2010, Sellers & Lotze 2015). \eqn{V=\mu\phi}{V=mu*phi}.}
##'      \item{genpois}{Generalized Poisson distribution (Consul & Famoye 1992). \eqn{V=\mu\exp(\eta)}{V=mu*exp(eta)}. (Note that Consul & Famoye (1992) define \eqn{\phi}{phi} differently.) Our implementation is taken from the \code{HMMpa} package, based on Joe and Zhu (2005) and implemented by Vitali Witowski.}
##'      \item{beta}{Beta distribution: parameterization of Ferrari and Cribari-Neto (2004)
##' and the \pkg{betareg} package (Cribari-Neto and Zeileis 2010); \eqn{V=\mu(1-\mu)/(\phi+1)}{V=mu*(1-mu)/(phi+1)}}
##'     \item{betabinomial}{Beta-binomial distribution: parameterized according to Morris (1997). \eqn{V=\mu(1-\mu)(n(\phi+n)/(\phi+1))}{V=mu*(1-mu)*(n*(phi+n)/(phi+1))}}
##'      \item{tweedie}{Tweedie distribution: \eqn{V=\phi\mu^p}{V=phi*mu^p}. The power parameter is restricted to the interval \eqn{1<p<2}. Code taken from the \code{tweedie} package, written by Peter Dunn.}
##' }
##' @references
##' \itemize{
##' \item Consul PC & Famoye F (1992). "Generalized Poisson regression model." Communications in Statistics: Theory and Methods 21:89–109.
##' \item Ferrari SLP, Cribari-Neto F (2004). "Beta Regression for Modelling Rates and Proportions." \emph{J. Appl. Stat.}  31(7), 799-815.
##' \item Hardin JW & Hilbe JM (2007). "Generalized linear models and extensions." Stata Press.
##' \item Huang A (2017). "Mean-parametrized Conway–Maxwell–Poisson regression models for dispersed counts." \emph{Statistical Modelling} 17(6), 1-22.
##' \item Joe H, Zhu R (2005). "Generalized Poisson Distribution: The Property of Mixture of Poisson and Comparison with Negative Binomial Distribution." \emph{Biometrical Journal} 47(2): 219–29. doi:10.1002/bimj.200410102.
##' \item Morris  W (1997). "Disentangling Effects of Induced Plant Defenses and Food Quantity on Herbivores by Fitting Nonlinear Models." \emph{American Naturalist} 150:299-327.
##' \item Sellers K & Lotze T (2015). "COMPoissonReg: Conway-Maxwell Poisson (COM-Poisson) Regression". R package version 0.3.5. https://CRAN.R-project.org/package=COMPoissonReg
##' \item Sellers K & Shmueli G (2010) "A Flexible Regression Model for Count Data." \emph{Annals of Applied Statistics} 4(2), 943–61. https://doi.org/10.1214/09-AOAS306.
##' \item Shonkwiler, J. S. (2016). "Variance of the truncated negative binomial distribution." \emph{Journal of Econometrics} 195(2), 209–210. doi:10.1016/j.jeconom.2016.09.002
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
#' @importFrom stats dnbinom pnbinom
#' @export	       	
truncated_nbinom2 <- function(link="log") {
    theta_errstr <- "theta (nbinom parameter) neither passed as an argument nor stored in enviroment"
    missing_theta <- "one" ## or "stop" or "na"
    r <- list(family="truncated_nbinom2",
              variance=function(mu,theta) {
                      if (missing(theta)) {
                          if (!exists(".Theta")) {
                              theta <- switch(missing_theta, one = 1, na = NA_real_, 
                                              stop = stop(theta_errstr))
                          }
                          else {
                              theta <- .Theta
                          }
                      }
                      a <- theta
                      c <- 0 ## truncation point
                      mu_star <- mu + (a+mu)*(c+1)*dnbinom(c+1,mu=mu,size=theta)/
                          (a*(1-pnbinom(c,mu=mu,size=theta)))
                      return(mu_star + c*(mu_star-mu) +mu_star*mu*(1+1/a)-mu_star^2)
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
                  if (exists("ziformula") && !ident(ziformula, ~0)) {
                      if (any(y < 0 | y >= 1)) {
                          stop("y values must be 0 <= y < 1")
                      }
                  } else {
                      if (any(y <= 0 | y >= 1)) 
                          stop("y values must be 0 < y < 1")
                  }
                  mustart <- y
              }))
    return(make_family(r,link))
}

## fixme: better name?

#' @rdname nbinom2
#' @export
## variance= (Wikipedia)
## n*alpha*beta*( alpha + beta + n )/ ((alpha+beta)^2*(alpha+beta+1))
## alpha = p*theta
## beta = (1-p)*theta
## -> n*p*(1-p)*theta^2*(theta+n)/(theta^2*(theta+1))
## =  n*p*(1-p)*(theta+n)/(theta+1)
##  *scaled* variance (dependence on mu only) is still just mu*(1-mu);
##  scaling is n*(theta+n)/(theta+1) (vs. simply n for the binomial)
betabinomial <- function(link="logit") {
    r <- list(family="betabinomial",
              variance=function(mu,phi) {
                  mu*(1-mu)
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

#' @export
#' @rdname nbinom2
ziGamma <- function(link="inverse") {
    g <- stats::Gamma(link=link)
    ## stats::Gamma does clever deparsing stuff ... need to work around it ...
    if (is.function(link)) {
        g$link <- deparse(substitute(link))
    } else g$link <- link
    ## modify initialization to allow zero values in zero-inflated cases
    g$initialize <- expression({
        if (exists("ziformula") && !ident(ziformula, ~0)) {
            if (any(y < 0)) stop("negative values not allowed for the 'Gamma' family with zero-inflation")
            } else {
                if (any(y <= 0)) stop("non-positive values not allowed for the 'Gamma' family")
            }
            n <- rep.int(1, nobs)
            mustart <- y
    })

   return(g)
}

