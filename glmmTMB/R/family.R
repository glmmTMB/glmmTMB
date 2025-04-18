## FIXME: I would like to use the following function instead of repeating
## the pattern, but I'm worried that lazy evaluation of arguments will
## cause all kinds of trouble
family_factory <- function(default_link, family, variance) {
    f <- function(link=default_link) {
        r <- list(family=family, link=link, variance=variance)
        r <- c(r,make.link(link))
        class(r) <- "family"
        return(r)
    }
    return(f)
}

## suppress code warnings for nbinom1/nbinom2; can't use .Theta <- NULL trick here ...
utils::globalVariables(".Theta")
utils::globalVariables(".Phi")

## see whether calling function has been called from glm.fit ...
in_glm_fit <- function() {
    up_two <- sys.calls()[[sys.nframe()-2]]
    identical(up_two[[1]], quote(glm.fit))
}

make_family <- function(x, link, needs_nonneg = FALSE, needs_int = FALSE) {
    if (is.character(link)) {
        x <- c(x, list(link=link), make.link(link))
    } else {
        x <- c(x, list(link=link$name), link)
    }
    ## stubs for Effect.default/glm.fit
    if (is.null(x$aic)) {
        x <- c(x,list(aic=function(...) NA_real_))
    }
    if (is.null(x$initialize)) {
        x <- c(x, list(initialize=
                           substitute(env = list(FAMILY=x$family,
                                                 needs_nonneg = needs_nonneg,
                                                 needs_int = needs_int),
            expr = expression({
            ## should handle log-links adequately
            mustart <- y+0.1
            if (needs_int) {
                if (any(abs(y - round(y)) > 0.001)) {
                    warning(gettextf("non-integer counts in a %s response variable", 
                                     FAMILY), domain = NA)
                }
            }
            if (needs_nonneg) {
                if (any(y < 0)) {
                    warning(gettextf("negative values in a %s response variable", 
                                     FAMILY), domain = NA)
                }
            }
            }))))
        ## strip one layer of protection
        ##   str == 'language expression ...' -> 'expression ...'
        x$initialize <- eval(x$initialize)

    }
        
    if (is.null(x$dev.resids)) {
        x <- c(x,list(dev.resids=function(y,mu,wt)  {
            if (in_glm_fit()) {
                ## can't return NA, glm.fit is unhappy
                ## glm() is called by the effects package
                ## cat("in glm.fit\n")
                return(rep(0,length(y)))
            } else {
                warning(warningCondition(
                    paste0("deviance residuals not defined ",
                          "for family ",
                          sQuote(x$family),
                          ": returning NA"),
                    class = c("dev_resids_undefined", "glmmTMB_warn")))
                return(rep(NA_real_, length(y)))
            }
        }))
    }

    class(x) <- "family"
    return(x)
}

## even better (?) would be to have a standalone list including
## name, default link, variance function, (optionally) initialize
## for each family


nbinom_errstr <- function(pname = "theta") {
    sprintf("%s (nbinom parameter) neither passed as an argument nor stored in enviroment", pname)
}
missing_disp <- "stop"

## find dispersion parameter in environment, if possible, or fall back
get_nbinom_disp <- function(disp, pname1 = ".Theta", pname2 = "theta") {
    if (is.null(disp)) {
        ## look in environment
        if (!exists(pname1, parent.frame())) {
            disp <- switch(missing_disp,
                           one=1,
                           na=NA_real_,
                           stop=stop(nbinom_errstr(pname2)))
        } else {
            disp <- get(pname1, parent.frame())
        }
    }
    assign(pname2, disp, parent.frame())
}

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
##' If specified, the dispersion model uses a log link; additional family parameters
##' (Student-t df, Tweedie power parameters, ordered beta cutpoints,
##' skew-normal skew parameters, etc.) use various link functions and
##' are accessible via \code{\link{family_params}}.
##' Denoting the variance as \eqn{V}, the dispersion parameter
##' as \eqn{\phi=\exp(\eta)}{phi=exp(eta)} (where \eqn{\eta}{eta} is the linear predictor from the dispersion model),
##' and the predicted mean as \eqn{\mu}{mu}:
##'  \describe{
##'      \item{gaussian}{(from base R): constant \eqn{V=\phi^2}{V=phi^2}}
##'      \item{Gamma}{(from base R) phi is the shape parameter. \eqn{V=\mu\phi}{V=mu*phi}}
##'       \item{ziGamma}{a modified version of \code{Gamma} that skips checks for zero values, allowing it to be used to fit hurdle-Gamma models}
##'      \item{nbinom2}{Negative binomial distribution: quadratic parameterization (Hardin & Hilbe 2007). \eqn{V=\mu(1+\mu/\phi) = \mu+\mu^2/\phi}{V=mu*(1+mu/phi) = mu+mu^2/phi}.}
##'      \item{nbinom1}{Negative binomial distribution: linear parameterization (Hardin & Hilbe 2007). \eqn{V=\mu(1+\phi)}{V=mu*(1+phi)}. \emph{Note} that the \eqn{phi} parameter has opposite meanings in the \code{nbinom1} and \code{nbinom2} families. In \code{nbinom1} overdispersion increases with increasing \code{phi} (the Poisson limit is \code{phi=0}); in \code{nbinom2} overdispersion decreases with increasing \code{phi} (the Poisson limit is reached as \code{phi} goes to infinity).}
##'      \item{nbinom12}{Negative binomial distribution: mixed linear/quadratic, as in the \code{DESeq2} package or as described by Lindén and Mäntyniemi (2011). \eqn{V=\mu(1+\phi+\mu/psi)}{V=mu*(1+phi+mu/psi)}. (In Lindén and Mäntyniemi's parameterization, \eqn{\omega = \phi}{omega=phi} and \eqn{\theta=1/\psi}{theta=1/psi}.) If a dispersion model is specified, it applies only to the linear (\code{phi}) term.}
##'      \item{truncated_nbinom2}{Zero-truncated version of nbinom2: variance expression from Shonkwiler 2016. Simulation code (for this and the other truncated count distributions) is taken from C. Geyer's functions in the \code{aster} package; the algorithms are described in \href{https://cran.r-project.org/package=aster/vignettes/trunc.pdf}{this vignette}.}
##'      \item{compois}{Conway-Maxwell Poisson distribution: parameterized with the exact mean (Huang 2017), which differs from the parameterization used in the \pkg{COMPoissonReg} package (Sellers & Shmueli 2010, Sellers & Lotze 2015). \eqn{V=\mu\phi}{V=mu*phi}.}
##'      \item{genpois}{Generalized Poisson distribution (Consul & Famoye 1992). \eqn{V=\mu\exp(\eta)}{V=mu*exp(eta)}. (Note that Consul & Famoye (1992) define \eqn{\phi}{phi} differently.) Our implementation is taken from the \code{HMMpa} package, based on Joe and Zhu (2005) and implemented by Vitali Witowski.}
##'      \item{beta}{Beta distribution: parameterization of Ferrari and Cribari-Neto (2004)
##' and the \pkg{betareg} package (Cribari-Neto and Zeileis 2010); \eqn{V=\mu(1-\mu)/(\phi+1)}{V=mu*(1-mu)/(phi+1)}}
##'     \item{betabinomial}{Beta-binomial distribution: parameterized according to Morris (1997). \eqn{V=\mu(1-\mu)(n(\phi+n)/(\phi+1))}{V=mu*(1-mu)*(n*(phi+n)/(phi+1))}}
##'      \item{tweedie}{Tweedie distribution: \eqn{V=\phi\mu^{power}}{V=phi*mu^power}. The power parameter is restricted to the interval \eqn{1<power<2}, i.e. the compound Poisson-gamma distribution. Code taken from the \code{tweedie} package, written by Peter Dunn. The power parameter (designated \code{psi} in the list of parameters) uses the link function \code{qlogis(psi-1.0)}; thus one can fix the power parameter to a specified value using \code{start = list(psi = qlogis(fixed_power-1.0)), map = list(psi = factor(NA))}.}
##'      \item{t_family}{Student-t distribution with adjustable scale and location parameters (also called a \href{https://en.wikipedia.org/wiki/Pearson_distribution#The_Pearson_type_VII_distribution}{Pearson type VII distribution}). The shape (degrees of freedom parameter) is fitted with a log link; it may be often be useful to fix the shape parameter using \code{start = list(psi = log(fixed_df)), map = list(psi = factor(NA))}.}
##'      \item{ordbeta}{Ordered beta regression from Kubinec (2022); fits continuous (e.g. proportion) data in the \emph{closed} interval [0,1]. Unlike the implementation in the \code{ordbeta} package, this family will not automatically scale the data. If your response variable is defined on the closed interval [a,b], transform it to [0,1] via \code{y_scaled <- (y-a)/(b-a)}.}
##'      \item{lognormal}{Log-normal, parameterized by the mean and standard deviation \emph{on the data scale}}
##'      \item{skewnormal}{Skew-normal, parameterized by the mean, standard deviation, and shape (Azzalini & Capitanio, 2014); constant \eqn{V=\phi^2}{V=phi^2}}
##' \item{bell}{Bell distribution (see Castellares et al 2018).
##' } 
##' }
##' @references
##' \itemize{
##' \item Azzalini A & Capitanio A (2014). "The skew-normal and related families." Cambridge: Cambridge University Press.
##' \item Castellares F, Ferrari SLP, & Lemonte AJ (2018) "On the Bell Distribution and Its Associated Regression Model for Count Data" Applied Mathematical Modelling 56: 172–85. \doi{10.1016/j.apm.2017.12.014}
##' \item Consul PC & Famoye F (1992). "Generalized Poisson regression model." Communications in Statistics: Theory and Methods 21:89–109.
##' \item Ferrari SLP, Cribari-Neto F (2004). "Beta Regression for Modelling Rates and Proportions." \emph{J. Appl. Stat.}  31(7), 799-815.
##' \item Hardin JW & Hilbe JM (2007). "Generalized linear models and extensions." Stata Press.
##' \item Huang A (2017). "Mean-parametrized Conway–Maxwell–Poisson regression models for dispersed counts." \emph{Statistical Modelling} 17(6), 1-22.
##' \item Joe H & Zhu R (2005). "Generalized Poisson Distribution: The Property of Mixture of Poisson and Comparison with Negative Binomial Distribution." \emph{Biometrical Journal} 47(2): 219–29. \doi{10.1002/bimj.200410102}.
##' \item Lindén, A & Mäntyniemi S. (2011). "Using the Negative Binomial Distribution to Model Overdispersion in Ecological Count Data." \emph{Ecology} 92 (7): 1414–21. \doi{10.1890/10-1831.1}.
##' \item Morris  W (1997). "Disentangling Effects of Induced Plant Defenses and Food Quantity on Herbivores by Fitting Nonlinear Models." \emph{American Naturalist} 150:299-327.
##' \item Kubinec R (2022). "Ordered Beta Regression: A Parsimonious, Well-Fitting Model for Continuous Data with Lower and Upper Bounds." \emph{Political Analysis}. doi:10.1017/pan.2022.20.
##' \item Sellers K & Lotze T (2015). "COMPoissonReg: Conway-Maxwell Poisson (COM-Poisson) Regression". R package version 0.3.5. https://CRAN.R-project.org/package=COMPoissonReg
##' \item Sellers K & Shmueli G (2010) "A Flexible Regression Model for Count Data." \emph{Annals of Applied Statistics} 4(2), 943–61. \doi{10.1214/09-AOAS306}.
##' \item Shonkwiler, J. S. (2016). "Variance of the truncated negative binomial distribution." \emph{Journal of Econometrics} 195(2), 209–210. \doi{10.1016/j.jeconom.2016.09.002}.
##' }
##' @export
##' @importFrom stats make.link
nbinom2 <- function(link="log") {
    r <- list(family="nbinom2",
              variance=function(mu, theta=NULL) {
                  get_nbinom_disp(theta, ".Theta", "theta")
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
              dev.resids = function (y, mu, wt, theta = NULL)  {
                  get_nbinom_disp(theta, ".Theta", "theta")
                  return(2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta))))
              })
    return(make_family(r,link))
}

#' @rdname nbinom2
#' @export
nbinom1 <- function(link="log") {
    r <- list(family="nbinom1",
              variance=function(mu, phi=NULL) {
                  get_nbinom_disp(phi, ".Phi", "phi")
                  return(mu*(1+phi))
              },
              initialize = expression({
                  if (any(y < 0))
                      stop("negative values not allowed for the negative binomial family")
                  n <- rep(1, nobs)
                  mustart <- y + (y == 0)/6
              }),
              dev.resids = function (y, mu, wt, phi = NULL)  {
                  get_nbinom_disp(phi, ".Phi", "phi")
                  ## convert phi to theta and use nbinom2 expression
                  ## V = mu*(1+phi) = mu*(1+mu/theta) -> theta = mu/phi
                  theta <- mu/phi
                  return(2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta))))
              })
    return(make_family(r,link))
}

#' @export
#' @rdname nbinom2
nbinom12 <- function(link="log") {
    r <- list(family="nbinom12",
              variance = function(mu, phi, psi) {
                  return(mu*(1+phi + mu/psi))
              }
              )
    return(make_family(r,link, needs_nonneg = TRUE, needs_int = TRUE))
}

#' @rdname nbinom2
#' @export
compois <- function(link="log") {
    r <- list(family="compois",
           variance=function(mu,phi) {
               if (length(phi)==1) phi <- rep(phi, length=length(mu))
               .Call("compois_calc_var", mu, 1/phi, PACKAGE="glmmTMB")
          })
    return(make_family(r, link, needs_nonneg = TRUE, needs_int = TRUE))
}

#' @rdname nbinom2
#' @export
truncated_compois <- function(link="log") {
    r <- list(family="truncated_compois",
           variance=function(mu,phi) {
             stop("variance for truncated compois family not yet implemented")
           })
    return(make_family(r,link, needs_nonneg = TRUE, needs_int = TRUE))
}

#' @rdname nbinom2
#' @export
genpois <- function(link="log") {
    r <- list(family="genpois",
           variance=function(mu,phi) {
               mu*phi
           })
    return(make_family(r,link, needs_nonneg = TRUE, needs_int = TRUE))
}

#' @rdname nbinom2
#' @export
truncated_genpois <- function(link="log") {
    r <- list(family="truncated_genpois",
           variance=function(mu,phi) {
             stop("variance for truncated genpois family not yet implemented")
          })
    return(make_family(r,link, needs_nonneg = TRUE, needs_int = TRUE))
}

#' @rdname nbinom2
#' @export
truncated_poisson <- function(link="log") {
	r <- list(family="truncated_poisson",
           variance=function(lambda) {
           (lambda+lambda^2)/(1-exp(-lambda)) - lambda^2/((1-exp(-lambda))^2)
           })
        return(make_family(r,link, needs_nonneg = TRUE, needs_int = TRUE))
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
    return(make_family(r,link, needs_nonneg = TRUE, needs_int = TRUE))
}

#' @rdname nbinom2
#' @export
truncated_nbinom1 <- function(link="log") {
    r <- list(family="truncated_nbinom1",
           variance=function(mu,alpha) {
               stop("variance for truncated nbinom1 family not yet implemented")
           })
    return(make_family(r,link, needs_nonneg = TRUE, needs_int = TRUE))
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
              })
              )
    return(make_family(r,link))
}

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
              variance = function(mu, phi) {
                  mu*(1-mu)
              },
              initialize = our_binom_initialize(binomial()$initialize))
    ## FIXME: should add needs_int = TRUE ??
    return(make_family(r,link))
}

#' @rdname nbinom2
#' @export
tweedie <- function(link="log") {
    r <- list(family="tweedie",
           variance = function(mu, phi, power) {
               phi * mu ^ power
         })
    return(make_family(r,link, needs_nonneg = TRUE))
}

#' @rdname nbinom2
#' @export
skewnormal <- function(link="identity") {
  r <- list(family="skewnormal",
            variance = function(phi) {
              phi^2
            })
  return(make_family(r,link))
}


#' @rdname nbinom2
#' @export
lognormal <- function(link="log") {
    r <- list(family="lognormal",
              variance=function(mu,phi) phi^2,
              initialize = expression({
                  if (exists("ziformula") && !ident(ziformula, ~0)) {
                      if (any(y < 0)) {
                          stop("y values must be >= 0")
                      }
                      mustart <- y + 0.1
                  } else {
                      if (any(y <= 0)) {
                          stop("y values must be > 0 (may be =0 if ziformula is specified)")
                      }
                  }
              })
              )
    return(make_family(r,link))
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

#' @export
#' @rdname nbinom2
t_family <- function(link="identity") {
    r <- list(family="t",
           variance=function(mu, phi) {
               rep(phi, length(mu))
           })
    return(make_family(r,link))
}


#' @export
#' @rdname nbinom2
ordbeta <- function(link="logit") {
    r <- list(family="ordbeta",
                            initialize=expression({
                                if (any(y < 0 | y > 1))
                                    stop("y values must be 0 <= y <= 1")
                                mustart <- y
                            }),
              ## from beta: not sure this is right ... ??
              variance=function(mu) { warning("ordbeta variance function untested"); mu*(1-mu) }
              )
    return(make_family(r,link))
}

#' @export
#' @rdname nbinom2
bell <- function(link="log") {
    ## can we get away with Suggests: gsl for this?
    ## or do we need Imports: ?
    if (!requireNamespace("gsl", quietly = TRUE)) {
        stop("the gsl package must be installed in order to use the Bell family")
    }
    r <- list(family="bell",
              variance = function(mu) {
                  mu*(1+gsl::lambert_W0(mu))
              }
              )
    ## Link <- list(linkfun = function(x) log(gsl::lambert_W0(x)),
    ##              linkinv = function(x) exp(x)*exp(exp(x)),
    ##              name = link ## mild hack to avoid make_family passing to make.link, which has hard-coded options
    ##              )
    return(make_family(r, link, needs_nonneg = TRUE, needs_int = TRUE))
}

