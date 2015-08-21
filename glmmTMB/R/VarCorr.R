

##' @S3method family glmmTMB
family.glmmTMB <- function(object, ...) {
    ## FIXME: this is wrong! -- family() must return a family, not a string !!!
    object$modelInfo$family
}

##' @S3method sigma glmmTMB
sigma.glmmTMB <- function(object, ...) {
    if(family(object) == "gaussian")
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

mkVC <- function(cor, sd, nms) {
    stopifnot(length(nms) == (nc <- length(cor)),  nc == length(sd),
              is.list(nms), is.list(cor), is.list(sd))

    do1cov <- function(sd, cor, n = length(sd)) sd * cor * rep(sd, each = n)

    if(FALSE)
        rr <- lapply(..., do1cov)
    ##

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


    structure(cov, stddev = sd, correlation = cor)


}

##' Extract variance and correlation components
##'
VarCorr.glmmTMB <- function(x, sigma = 1, rdig = 3)# <- 3 args from nlme
{

  ## FIXME:: add type=c("varcov","sdcorr","logs" ?)
    stopifnot(is.numeric(sigma), length(sigma) == 1)

    xrep <- x$obj$env$report()

    reT <- x$modelInfo$reTerms
    ## $ reTrms :List of 3
    ##   ..$ condList:List of 2
    ##   .. ..$ cnms :List of 1
    ##   .. .. ..$ Subject: chr [1:2] "(Intercept)" "age"

    ## ans <- structure(val, stddev = stddev, correlation = corr)
    vc.cond <- mkVC(cor = xrep$corr,  sd = xrep$sd,   cnms = reT$condList$cnms)
    vs.zi   <- mkVC(cor = xrep$corzi, sd = xrep$sdzi, cnms = reT$ ziList $cnms)

    ## unfinished
    structure(ans, sc = sc)


   structure(mkVarCorr(sigma, cnms = cnms, nc = nc, theta = x@theta,
			nms = { fl <- x@flist; names(fl)[attr(fl, "assign")]}),
	      useSc = as.logical(x@devcomp$dims[["useSc"]]),
	      class = "VarCorr.merMod")
}


##' @S3method print VarCorr.merMod
print.VarCorr.glmmTMB <- function(x, digits = max(3, getOption("digits") - 2),
		   comp = "Std.Dev.", formatter = format, ...) {
    print(formatVC(x, digits = digits, comp = comp, formatter = formatter), quote = FALSE, ...)
    invisible(x)
}

