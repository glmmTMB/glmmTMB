##' Extract the modes of the random effects
##'
##' A generic function to extract the conditional modes of the random effects
##' from a fitted model object.  For linear mixed models the conditional modes
##' of the random effects are also the conditional means.
##'
##' If grouping factor i has k levels and j random effects per level the ith
##' component of the list returned by \code{ranef} is a data frame with k rows
##' and j columns.  If \code{condVar} is \code{TRUE} the \code{"postVar"}
##' attribute is an array of dimension j by j by k.  The kth face of this array
##' is a positive definite symmetric j by j matrix.  If there is only one
##' grouping factor in the model the variance-covariance matrix for the entire
##' random effects vector, conditional on the estimates of the model parameters
##' and on the data will be block diagonal and this j by j matrix is the kth
##' diagonal block.  With multiple grouping factors the faces of the
##' \code{"postVar"} attributes are still the diagonal blocks of this
##' conditional variance-covariance matrix but the matrix itself is no longer
##' block diagonal.
##' @name ranef
##' @aliases ranef ranef.merMod
##' @param object an object of a class of fitted models with random effects,
##' typically an \code{"\linkS4class{merMod}"} object.
##' @param condVar an optional logical argument indicating if the conditional
##' variance-covariance matrices of the random effects should be added as an attribute.
##' @param postVar a (deprecated) synonym for \code{condVar}
##' @param drop an optional logical argument indicating components of the return
##' value that would be data frames with a single column, usually a column
##' called \sQuote{\code{(Intercept)}}, should be returned as named vectors.
##' @param whichel an optional character vector of names of grouping factors for
##' which the random effects should be returned.  Defaults to all the grouping
##' factors.
##' @param \dots some methods for this generic function require additional
##' arguments.
##' @return A list of data frames, one for each grouping factor for the random
##' effects.  The number of rows in the data frame is the number of levels of
##' the grouping factor.  The number of columns is the dimension of the random
##' effect associated with each level of the factor.
##'
##' If \code{condVar} is \code{TRUE} each of the data frames has an attribute
##' called \code{"postVar"} which is a three-dimensional array with symmetric
##' faces.
##'
##' When \code{drop} is \code{TRUE} any components that would be data frames of
##' a single column are converted to named numeric vectors.
##' @note To produce a \dQuote{caterpillar plot} of the random effects apply
##' \code{\link[lattice:xyplot]{dotplot}} to the result of a call to
##' \code{ranef} with \code{condVar = TRUE}.
##' @examples
##' fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
##' fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
##' fm3 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
##' ranef(fm1)
##' str(rr1 <- ranef(fm1, condVar = TRUE))
##' dotplot(rr1)  ## default
##' ## specify free scales in order to make Day effects more visible
##' dotplot(rr1,scales = list(x = list(relation = 'free')))[["Subject"]]
##' if(FALSE) { ##-- condVar=TRUE is not yet implemented for multiple terms -- FIXME
##' str(ranef(fm2, condVar = TRUE))
##' }
##' op <- options(digits = 4)
##' ranef(fm3, drop = TRUE)
##' options(op)
##' @keywords models methods
##' @method ranef merMod
##' @export
ranef.merMod <- function(object, condVar = FALSE, drop = FALSE,
			 whichel = names(ans), postVar = FALSE, ...)
{
    if (!missing(postVar) && missing(condVar)) {
        warning(sQuote("postVar")," is deprecated: please use ",
                sQuote("condVar")," instead")
        condVar <- postVar
    }
    ans <- c(matrix(unlist(getME(object,"b")))) ## == object@pp$b(1.)
    if (!is.null(object@flist)) {
	## evaluate the list of matrices
	levs <- lapply(fl <- object@flist, levels)
	asgn <- attr(fl, "assign")
	cnms <- object@cnms
	nc <- vapply(cnms, length, 1L)
	nb <- nc * vapply(levs, length, 1L)[asgn]
	nbseq <- rep.int(seq_along(nb), nb)
	ml <- split(ans, nbseq)
	for (i in seq_along(ml))
	    ml[[i]] <- matrix(ml[[i]], ncol = nc[i], byrow = TRUE,
			      dimnames = list(NULL, cnms[[i]]))
	## create a list of data frames corresponding to factors
	ans <- lapply(seq_along(fl),
		      function(i)
		      data.frame(do.call(cbind, ml[asgn == i]),
				 row.names = levs[[i]],
				 check.names = FALSE))
	names(ans) <- names(fl)
					# process whichel
	stopifnot(is(whichel, "character"))
	whchL <- names(ans) %in% whichel
	ans <- ans[whchL]

	if (condVar) {
            sigsqr <- sigma(object)^2
            rp <- rePos$new(object)
            if(any(sapply(rp$terms, length) > 1)){
                # TODO: actually use condVar and then convert back to array format
                warning("conditional variances not currently available via ",
                        "ranef when there are multiple terms per factor")
            } else{
                vv <- .Call(merPredDcondVar, object@pp$ptr(), as.environment(rp))
                for (i in names(ans)) ## seq_along(ans))
                    attr(ans[[i]], "postVar") <- vv[[i]] * sigsqr
            }
	}
	if (drop)
	    ans <- lapply(ans, function(el)
		      {
			  if (ncol(el) > 1) return(el)
			  pv <- drop(attr(el, "postVar"))
			  el <- drop(as.matrix(el))
			  if (!is.null(pv))
                              attr(el, "postVar") <- pv
			  el
		      })
     	class(ans) <- "ranef.mer"
    }
    ans
}## ranef.merMod

print.ranef.mer <- function(x, ...) {
    print(unclass(x), ...)
    if(any(has.pv <- vapply(x, function(el)
			    !is.null(attr(el, "postVar")), NA)))
	cat('with conditional variances for',
	    paste(dQuote(names(x)[has.pv]), sep=", "), "\n")
    invisible(x)
}
