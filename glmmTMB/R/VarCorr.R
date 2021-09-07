## returns a true family() object iff one was given
## to glmmTMB() in the first place ....
##' @importFrom stats family
##' @export
family.glmmTMB <- function(object, ...) {
    object$modelInfo$family
}

## don't quite know why this (rather than just ...$parList()) is
## necessary -- used in ranef.glmmTMB and sigma.glmmTMB
getParList <- function(object) {
    object$obj$env$parList(object$fit$par, object$fit$parfull)
}


##' Extract residual standard deviation or dispersion parameter
##'
##' For Gaussian models, \code{sigma} returns the value of the residual
##' standard deviation; for other families, it returns the
##' dispersion parameter, \emph{however it is defined for that
##' particular family}. See details for each family below.
##'
##' @details
##'  The value returned varies by family:
##'  \describe{
##'      \item{gaussian}{returns the \emph{maximum likelihood} estimate
##'          of the standard deviation (i.e., smaller than the results of
##'                                  \code{sigma(lm(...))} by a factor of (n-1)/n)}
##'      \item{nbinom1}{returns a dispersion parameter
##'          (usually denoted \eqn{\alpha}{alpha} as in Hardin and Hilbe (2007)):
##'          such that the variance equals \eqn{\mu(1+\alpha)}{mu(1+alpha)}.}
##'      \item{nbinom2}{returns a dispersion parameter
##'          (usually denoted \eqn{\theta}{theta} or \eqn{k}); in contrast to
##'          most other families, larger \eqn{\theta}{theta} corresponds to a \emph{lower}
##'          variance which is \eqn{\mu(1+\mu/\theta)}{mu(1+mu/theta)}.}
##'      \item{Gamma}{Internally, glmmTMB fits Gamma responses by fitting a mean
##'          and a shape parameter; sigma is estimated as (1/sqrt(shape)),
##'          which will typically be close (but not identical to) that estimated
##'          by \code{stats:::sigma.default}, which uses sqrt(deviance/df.residual)}
##'      \item{beta}{returns the value of \eqn{\phi}{phi},
##'          where the conditional variance is \eqn{\mu(1-\mu)/(1+\phi)}{mu*(1-mu)/(1+phi)}
##'          (i.e., increasing \eqn{\phi}{phi} decreases the variance.)
##'          This parameterization follows Ferrari and Cribari-Neto (2004)
##'          (and the \code{betareg} package):}
##'      \item{betabinomial}{This family uses the same parameterization (governing
##'           the Beta distribution that underlies the binomial probabilities) as \code{beta}.}
##'      \item{genpois}{returns the index of dispersion \eqn{\phi^2}{phi^2},
##'           where the variance is \eqn{\mu\phi^2}{mu*phi^2} (Consul & Famoye 1992)}
##'      \item{compois}{returns the value of \eqn{1/\nu}{1/nu},
##'           When \eqn{\nu=1}{nu=1}, compois is equivalent to the Poisson distribution.
##'           There is no closed form equation for the variance, but
##'           it is approximately undersidpersed when \eqn{1/\nu <1}{1/nu <1}
##'           and approximately oversidpersed when \eqn{1/\nu >1}{1/nu>1}.
##'           In this implementation, \eqn{\mu}{mu} is exactly the mean (Huang 2017), which
##'           differs from the COMPoissonReg package (Sellers & Lotze 2015).}
##'      \item{tweedie}{returns the value of \eqn{\phi}{phi},
##'           where the variance is \eqn{\phi\mu^p}{phi*mu^p}.
##'           The value of \eqn{p} can be extracted using the internal
##'           function \code{glmmTMB:::.tweedie_power}.}
##' }
##'
##'  The most commonly used GLM families
##'  (\code{binomial}, \code{poisson}) have fixed dispersion parameters which are
##'  internally ignored.
##'
##' @references
##' \itemize{
##' \item Consul PC, and Famoye F (1992). "Generalized Poisson regression model. Communications in Statistics: Theory and Methods" 21:89–109.
##' \item Ferrari SLP, Cribari-Neto F (2004). "Beta Regression for Modelling Rates and Proportions." \emph{J. Appl. Stat.}  31(7), 799-815.
##' \item Hardin JW & Hilbe JM (2007). "Generalized linear models and extensions." Stata press.
##' \item Huang A (2017). "Mean-parametrized Conway–Maxwell–Poisson regression models for dispersed counts. " \emph{Statistical Modelling} 17(6), 1-22.
##' \item Sellers K & Lotze T (2015). "COMPoissonReg: Conway-Maxwell Poisson (COM-Poisson) Regression". R package version 0.3.5. https://CRAN.R-project.org/package=COMPoissonReg
##' }
##' @aliases sigma
##' @param object a \dQuote{glmmTMB} fitted object
##' @param \dots (ignored; for method compatibility)
## Import generic and re-export
## note the following line is hacked in Makefile/namespace-update to ...
## if(getRversion()>='3.3.0') importFrom(stats, sigma) else importFrom(lme4,sigm
## also see <https://github.com/klutometis/roxygen/issues/371>
##' @rawNamespace if(getRversion()>='3.3.0') importFrom(stats, sigma) else importFrom(lme4,sigma)
##  n.b. REQUIRES roxygen2 >= 5.0
## @importFrom lme4 sigma
##' @export sigma
##' @method sigma glmmTMB
##' @export
sigma.glmmTMB <- function(object, ...) {
    pl <- getParList(object)
    ff <- object$modelInfo$family$family
    if (!usesDispersion(ff)) return(1.)
    if (length(pl$betad)>1) return(NA)
    switch(family(object)$family,
           gaussian=exp(0.5*pl$betad),
           Gamma=exp(-0.5*pl$betad),
           exp(pl$betad))
}


mkVC <- function(cor, sd, cnms, sc, useSc) {
    stopifnot(length(cnms) == (nc <- length(cor)),  nc == length(sd),
              is.list(cnms), is.list(cor), is.list(sd),
              is.character(nnms <- names(cnms)), nzchar(nnms))
    ##
    ## FIXME: do we want this?  Maybe not.
    ## Potential problem: the names of the elements of the VarCorr() list
    ##  are not necessarily unique (e.g. fm2 from example("lmer") has *two*
    ##  Subject terms, so the names are "Subject", "Subject".  The print method
    ##  for VarCorrs handles this just fine, but it's a little awkward if we
    ##  want to dig out elements of the VarCorr list ... ???
    if (anyDuplicated(nnms))
        nnms <- make.names(nnms, unique = TRUE)
    ##
    ## cov :=  F(sd, cor) :
    do1cov <- function(sd, cor, n = length(sd)) {
        sd * cor * rep(sd, each = n)
    }
    docov <- function(sd,cor,nm) {
        ## diagonal model:
        diagmodel <- identical(dim(cor),c(0L,0L))
        if (diagmodel) cor <- diag(length(sd))
        cov <- do1cov(sd, cor)
        names(sd) <- nm
        dimnames(cov) <- dimnames(cor) <- list(nm,nm)
        structure(cov,stddev=sd,correlation=cor)
    }
    ss <- setNames(mapply(docov,sd,cor,cnms,SIMPLIFY=FALSE),nnms)
    ## ONLY first element -- otherwise breaks formatVC
    ## FIXME: do we want a message/warning here, or elsewhere,
    ##   when the 'Residual' var parameters are truncated?
    attr(ss,"sc") <- sc[1]
    attr(ss,"useSc") <- useSc
    ss
}

##' Extract variance and correlation components
##'
##' @aliases VarCorr
##' @param x a fitted \code{glmmTMB} model
##' @param sigma residual standard deviation (usually set automatically from internal information)
##' @param extra arguments (for consistency with generic method)
##' @importFrom nlme VarCorr
## and re-export the generic:
##' @export VarCorr
##' @export
##' @examples
##' ## Comparing variance-covariance matrix with manual computation
##' data("sleepstudy",package="lme4")
##' fm4 <- glmmTMB(Reaction ~ Days + (Days|Subject), sleepstudy)
##' VarCorr(fm4)[[c("cond","Subject")]]
##' ## hand calculation
##' pars <- getME(fm4,"theta")
##' ## construct cholesky factor
##' L <- diag(2)
##' L[lower.tri(L)] <- pars[-(1:2)]
##' C <- crossprod(L)
##' diag(C) <- 1
##' sdvec <- exp(pars[1:2])
##' (V <- outer(sdvec,sdvec) * C)
##' @details For an unstructured variance-covariance matrix, the internal parameters
##' are structured as follows: the first n parameters are the log-standard-deviations,
##' while the remaining n(n-1)/2 parameters are the elements of the Cholesky factor
##' of the correlation matrix, filled in column-wise order
##' (see the \href{http://kaskr.github.io/adcomp/classUNSTRUCTURED__CORR__t.html}{TMB documentation}
##' for further details).
##' @keywords internal
VarCorr.glmmTMB <- function(x, sigma = 1, ... )
{
    ## FIXME:: add type=c("varcov","sdcorr","logs" ?)
    ## FIXME:: do we need 'sigma' any more (now that nlme generic
    ##         doesn't have it?)
    stopifnot(is.numeric(sigma), length(sigma) == 1)
    xrep <- x$obj$env$report(x$fit$parfull)
    reT <- x$modelInfo$reTrms
    reS <- x$modelInfo$reStruc
    familyStr <- family(x)$family
    useSc <- if (missing(sigma)) {
        ## *only* report residual variance for Gaussian family ...
        ## *not* usesDispersion(familyStr)
        sigma <- sigma(x)
        familyStr=="gaussian" && !zeroDisp(x)
    } else TRUE
    vc.cond <- vc.zi <- NULL
    if(length(cn <- reT$cond$cnms)) {
        vc.cond <- mkVC(cor = xrep$corr,  sd = xrep$sd,   cnms = cn,
                        sc = sigma, useSc = useSc)
        for (i in seq_along(vc.cond)) {
            attr(vc.cond[[i]],"blockCode") <- reS$condReStruc[[i]]$blockCode
        }
    }
    if(length(cn <- reT$zi$cnms)) {
        vc.zi <- mkVC(cor = xrep$corrzi, sd = xrep$sdzi, cnms = cn,
                        sc = sigma, useSc = useSc)
        for (i in seq_along(vc.zi)) {
            attr(vc.zi,"blockCode") <- reS$ziReStruc[[i]]$blockCode
        }
    }
    structure(list(cond = vc.cond, zi = vc.zi),
	      sc = usesDispersion(familyStr), ## 'useScale'
	      class = "VarCorr.glmmTMB")
}

##' Printing The Variance and Correlation Parameters of a \code{glmmTMB}
##' @method print VarCorr.glmmTMB
##' @export
## don't importFrom lme4 formatVC; use our own formatV instead!
##  document as it is a method with "surprising arguments":
##' @param x a result of \code{\link{VarCorr}(<glmmTMB>)}.
##' @param digits number of significant digits to use.
##' @param comp a string specifying the component to format and print.
##' @param formatter a \code{\link{function}}.
##' @param ... optional further arguments, passed the next \code{\link{print}} method.
print.VarCorr.glmmTMB <- function(x, digits = max(3, getOption("digits") - 2),
				  comp = "Std.Dev.", formatter = format, ...)
{
    for (cc in names(x))  if(!is.null(x[[cc]])) {
	cat(sprintf("\n%s:\n", cNames[[cc]]))
	print(formatVC(x[[cc]],
		       digits = digits, comp = comp, formatter = formatter),
	      quote = FALSE, ...)
    }
    invisible(x)
}

## original from lme4
## had to be extended/modified to deal with glmmTMB special cases
##' "format()" the 'VarCorr' matrix of the random effects -- for
##' print()ing and show()ing
##'
##' @title Format the 'VarCorr' Matrix of Random Effects
##' @param varcor a \code{\link{VarCorr}} (-like) matrix with attributes.
##' @param digits the number of significant digits.
##' @param comp character vector of length one or two indicating which
##' columns out of "Variance" and "Std.Dev." should be shown in the
##' formatted output.
##' @param formatter the \code{\link{function}} to be used for
##' formatting the standard deviations and or variances (but
##' \emph{not} the correlations which (currently) are always formatted
##' as "0.nnn"
##' @param useScale whether to report a scale parameter (e.g. residual standard deviation)
##' @param ... optional arguments for \code{formatter(*)} in addition
##' to the first (numeric vector) and \code{digits}.
##' @return a character matrix of formatted VarCorr entries from \code{varc}.
##' @importFrom methods as
formatVC <- function(varcor, digits = max(3, getOption("digits") - 2),
		     comp = "Std.Dev.", formatter = format,
                     useScale = attr(varcor, "useSc"),
                     ...)
{
    c.nms <- c("Groups", "Name", "Variance", "Std.Dev.")
    avail.c <- c.nms[-(1:2)]
    if(anyNA(mcc <- pmatch(comp, avail.c)))
	stop("Illegal 'comp': ", comp[is.na(mcc)])
    nc <- length(colnms <- c(c.nms[1:2], (use.c <- avail.c[mcc])))
    if(length(use.c) == 0)
	stop("Must show either standard deviations or variances")
    formatCor <- function(x,maxlen=0) {
        ## x: correlation matrix
        ## maxlen: max number of RE std devs per term;
        ##         really a *minimum* length here! pad to (maxlen)
        ##         columns as necessary
        x <- as(x, "matrix")
        dig <- max(2, digits - 2) # use 'digits' !
        ## n.b. not using formatter() for correlations
        cc <- format(round(x, dig), nsmall = dig)
        cc[!lower.tri(cc)] <- ""  ## empty lower triangle
        nr <- nrow(cc)
        if (nr < maxlen) {
            cc <- cbind(cc, matrix("", nr, maxlen-nr))
        }
        return(cc)
    }
    getCovstruct <- function(x) {
        n <- names(.valid_covstruct)[match(attr(x,"blockCode"),
                                           .valid_covstruct)]
        if (length(n)==0) n <- "us"  ## unstructured v-cov (default)
        return(n)
    }
    getCorSD <- function(x,type="stddev",maxlen=0) {
        r <- attr(x,type) ## extract stddev *or* correlation from x
        if (type=="correlation") {
            ## transform to char matrix
            r <- formatCor(r,maxlen)
            ## drop last column (will be blank since we blanked out
            ##  the upper triangle + diagonal)
            if (ncol(r)>maxlen)
                r <- r[, -ncol(r), drop = FALSE]
        }
        covstruct <- getCovstruct(x)
        if (covstruct %in% c("ar1","cs")) {
            r <- switch(type,
                        stddev=r[1],
                        ## select lag-1 correlation
                        ## upper tri has been erased in formatCor() ...
                        correlation=paste(r[2,1],sprintf("(%s)",covstruct))
                        )
        }
        return(r)
    }  ## getCorSD

    ## get std devs:
    reStdDev <- lapply(varcor, getCorSD)
    ## need correlations if
    useCor <- (sapply(varcor,getCovstruct)!="us" |
               sapply(reStdDev,length)>1)
    cnms <- Map(function(x,n) colnames(x)[seq(n)], varcor, lengths(reStdDev))

    if(useScale) {
        reStdDev <- c(reStdDev,
                      list(Residual = unname(attr(varcor, "sc"))))
     }

    reLens <- lengths(reStdDev)
    nr <- sum(reLens)
    reMat <- array('', c(nr, nc), list(rep.int('', nr), colnms))
    reMat[1+cumsum(reLens)-reLens, "Groups"] <- names(reLens)

    reMat[,"Name"] <- c(unlist(cnms), if(useScale) "")
    if("Variance" %in% use.c)
	reMat[,"Variance"] <- formatter(unlist(reStdDev)^2, digits = digits, ...)
    if("Std.Dev." %in% use.c)
	reMat[,"Std.Dev."] <- formatter(unlist(reStdDev),   digits = digits, ...)
    if (any(useCor)) {
        ## get corrs
	maxlen <- max(reLens)
	corr <-
	    do.call(rbind,lapply(varcor, getCorSD,
                                         type="correlation", maxlen=maxlen))
        ## add blank values as necessary
	if (nrow(corr) < nrow(reMat))
	    corr <- rbind(corr, matrix("", nrow(reMat) - nrow(corr), ncol(corr)))
	colnames(corr) <- c("Corr", rep.int("", max(0L, ncol(corr)-1L)))
	cbind(reMat, corr)
    } else reMat
}
