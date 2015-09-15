##' returns a true family() object iff one was given
##' to glmmTMB() in the first place ....
##' @importFrom stats family
##' @export family
##' @keywords internal
family.glmmTMB <- function(object, ...) {
    object$modelInfo$family
}

## Import generic and re-export
## note the following line is hacked in Makefile/namespace-update to ...
## if(getRversion()>='3.3.0') importFrom(stats, sigma) else importFrom(lme4,sigma)
## also see <https://github.com/klutometis/roxygen/issues/371>
##' @importFrom lme4 sigma
##' @export sigma

##' @method sigma glmmTMB
##' @export
##' @keywords internal
sigma.glmmTMB <- function(object, ...) {
    if(family(object)$family == "gaussian")
        exp( .5 * object$obj$env$parList()$betad ) # betad is  log(sigma ^ 2)
    else 1.
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
        ## FIXME: what should be in cor for a 1x1 diag model?
        if (identical(dim(cor),c(0L,0L))) cor <- matrix(1)
        cov <- do1cov(sd, cor)
        names(sd) <- nm
        dimnames(cov) <- dimnames(cor) <- list(nm,nm)
        structure(cov,stddev=sd,correlation=cor)
    }
    ss <- setNames(mapply(docov,sd,cor,cnms),nnms)
    attr(ss,"sc") <- sc
    attr(ss,"useSc") <- useSc
    ss
}


##' Extract variance and correlation components
##'
##' @importFrom nlme VarCorr
## and re-export the generic:
##' @export VarCorr
##' @export
VarCorr.glmmTMB <- function(x, sigma = 1, rdig = 3)# <- 3 args from nlme
{
    ## FIXME:: add type=c("varcov","sdcorr","logs" ?)
    stopifnot(is.numeric(sigma), length(sigma) == 1)
    xrep <- x$obj$env$report()
    reT <- x$modelInfo$reTrms
    useSc <- if (missing(sigma)) {
        sigma <- sigma(x)
        usesDispersion(family(x)$family)
    } else TRUE

    vc.cond <- if(length(cn <- reT$condList$cnms))
        mkVC(cor = xrep$corr,  sd = xrep$sd,   cnms = cn, sc = sigma, useSc = useSc)
    vc.zi   <- if(length(cn <- reT$  ziList$cnms))
        mkVC(cor = xrep$corzi, sd = xrep$sdzi, cnms = cn)
    structure(list(cond = vc.cond, zi = vc.zi),
	      sc = usesDispersion(x), ## 'useScale'
	      class = "VarCorr.glmmTMB")
}

##' @title Printing The Variance and Correlation Parameters of a \code{glmmTMB}
##' @method print VarCorr.glmmTMB
##' @export
##' @importFrom lme4 formatVC
##              ^^^^ github version >= 2015-09-05
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

