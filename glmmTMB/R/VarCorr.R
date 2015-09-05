##' returns a true family() object iff one was given
##' to glmmTMB() in the first place ....
##' @importFrom stats family
##' @export
##' @keywords internal
family.glmmTMB <- function(object, ...) {
    object$modelInfo$family
}

## note the following line is hacked in Makefile/namespace-update to ...
## if(getRversion()>='3.3.0') importFrom(stats, sigma) else importFrom(lme4,sigma)
## also see <https://github.com/klutometis/roxygen/issues/371>
##' @importFrom lme4 sigma
##' @export sigma
##' @export
##' @keywords internal
sigma.glmmTMB <- function(object, ...) {
    if(family(object)$family == "gaussian")
        exp( .5 * object$obj$env$parList()$betad ) # betad is  log(sigma ^ 2)
    else 1.
}


##' Make variance and correlation matrices from \code{theta}
##'
##' @param sc scale factor (residual standard deviation)
##' @param cnms component names
##' @param nc numeric vector: number of terms in each RE component
##' @param theta theta vector (lower-triangle of Cholesky factors)
##' @param nms component names (FIXME: nms/cnms redundant: nms=names(cnms)?)
##' @seealso \code{\link{VarCorr}}
##' @return A matrix
##' @export
##' @keywords internal
mkVarCorr <- function(sc, cnms, nc, theta, nms) {
    ncseq <- seq_along(nc)
    thl <- split(theta, rep.int(ncseq, (nc * (nc + 1))/2))
    if(!all(nms == names(cnms))) ## the above FIXME
	warning("nms != names(cnms)  -- whereas lme4-authors thought they were --\n",
		"Please report!", immediate. = TRUE)
    ans <- lapply(ncseq, function(i)
	      {
		  ## Li := \Lambda_i, the i-th block diagonal of \Lambda(\theta)
		  Li <- diag(nrow = nc[i])
		  Li[lower.tri(Li, diag = TRUE)] <- thl[[i]]
		  rownames(Li) <- cnms[[i]]
		  ## val := \Sigma_i = \sigma^2 \Lambda_i \Lambda_i', the
		  val <- tcrossprod(sc * Li) # variance-covariance
		  stddev <- sqrt(diag(val))
		  corr <- t(val / stddev)/stddev
		  diag(corr) <- 1
		  structure(val, stddev = stddev, correlation = corr)
	      })
    if(is.character(nms)) {
	## FIXME: do we want this?  Maybe not.
	## Potential problem: the names of the elements of the VarCorr() list
	##  are not necessarily unique (e.g. fm2 from example("lmer") has *two*
	##  Subject terms, so the names are "Subject", "Subject".  The print method
	##  for VarCorrs handles this just fine, but it's a little awkward if we
	##  want to dig out elements of the VarCorr list ... ???
	if (anyDuplicated(nms))
	    nms <- make.names(nms, unique = TRUE)
	names(ans) <- nms
    }
    structure(ans, sc = sc)
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
    do1cov <- function(sd, cor, n = length(sd)) sd * cor * rep(sd, each = n)
    ##
    ss <- setNames(lapply(seq_len(nc), function(i) {
        cov <- do1cov(sd = (sd.i <- sd[[i]]), cor = (cor.i <- cor[[i]]))
        names(sd.i) <- nm.i <- cnms[[i]]
        dimnames(cov) <- dimnames(cor.i) <- list(nm.i, nm.i)
        structure(cov, stddev = sd.i, correlation = cor.i)
    }), nnms)
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
    if (missing(sigma)) {
        sigma <- sigma(x)
        useSc <- usesDispersion(family(x)$family)
    } else {
        useSc <- TRUE
    }
    vc.cond <- if(length(cn <- reT$condList$cnms)) {
        mkVC(cor = xrep$corr,  sd = xrep$sd,   cnms = cn,
             sc = sigma, useSc = useSc)
        
    }
    vc.zi   <- if(length(cn <- reT$  ziList$cnms))
                   mkVC(cor = xrep$corzi, sd = xrep$sdzi, cnms = cn)
    structure(list(cond = vc.cond, zi = vc.zi),
	      class = "VarCorr.glmmTMB")
}

##' @title Printing The Variance and Correlation Parameters of a \code{glmmTMB}
##' @method print VarCorr.glmmTMB
##' @export
##  document as it is a method with "surprising arguments":
##' @param x a result of \code{\link{VarCorr}(<glmmTMB>)}.
##' @param digits number of significant digits to use
##' @param comp a string specifying the component to format and print
##' @param formatter a \code{\link{function}}
##' @param ...
print.VarCorr.glmmTMB <- function(x, digits = max(3, getOption("digits") - 2),
		   comp = "Std.Dev.", formatter = format, ...) {
    for (cc in names(x)) {
        if (!is.null(term <- x[[cc]])) {
            cat("\n",cNames[[cc]],"\n",sep="")
            print(formatVC(term,
                           digits = digits, comp = comp, formatter = formatter), quote = FALSE, ...)
        }
    }
    invisible(x)
}

