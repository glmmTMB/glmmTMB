## backward compat (copied from lme4)
if((Rv <- getRversion()) < "3.2.1") {
    lengths <- function (x, use.names = TRUE) vapply(x, length, 1L, USE.NAMES = use.names)
}
rm(Rv)

## generate a list with names equal to values
## See also: \code{tibble::lst}, \code{Hmisc::llist}
namedList <- function (...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L)))
        nm <- snm
    if (any(nonames <- nm == ""))
        nm[nonames] <- snm[nonames]
    setNames(L, nm)
}

## Sparse Schur complement (Marginal of precision matrix)
##' @importFrom Matrix Cholesky solve
GMRFmarginal <- function(Q, i, ...) {
    ind <- seq_len(nrow(Q))
    i1 <- (ind)[i]
    i0 <- setdiff(ind, i1)
    if (length(i0) == 0)
        return(Q)
    Q0 <- as(Q[i0, i0, drop = FALSE], "symmetricMatrix")
    L0 <- Cholesky(Q0, ...)
    ans <- Q[i1, i1, drop = FALSE] - Q[i1, i0, drop = FALSE] %*%
        solve(Q0, as.matrix(Q[i0, i1, drop = FALSE]))
    ans
}

parallel_default <- function(parallel=c("no","multicore","snow"),ncpus=1) {
    ##  boilerplate parallel-handling stuff, copied from lme4
    if (missing(parallel)) parallel <- getOption("profile.parallel", "no")
    parallel <- match.arg(parallel)
    do_parallel <- (parallel != "no" && ncpus > 1L)
    if (do_parallel && parallel == "multicore" &&
        .Platform$OS.type == "windows") {
        warning("no multicore on Windows, falling back to non-parallel")
        parallel <- "no"
    }
    return(list(parallel=parallel,do_parallel=do_parallel))
}

##' translate vector of correlation parameters to correlation values
##' @param theta vector of internal correlation parameters
##' @return a vector of correlation values
##' @details This function follows the definition at \url{http://kaskr.github.io/adcomp/classUNSTRUCTURED__CORR__t.html}:
##' if \eqn{L} is the lower-triangular matrix with 1 on the diagonal and the correlation parameters in the lower triangle, then the correlation matrix is defined as \eqn{\Sigma = D^{-1/2} L L^\top D^{-1/2}}{Sigma = sqrt(D) L L' sqrt(D)}, where \eqn{D = \textrm{diag}(L L^\top)}{D = diag(L L')}. For a single correlation parameter \eqn{\theta_0}{theta0}, this works out to \eqn{\rho = \theta_0/\sqrt{1+\theta_0^2}}{rho = theta0/sqrt(1+theta0^2)}. The function returns the elements of the lower triangle of the correlation matrix, in column-major order.
##' @examples
##' th0 <- 0.5
##' stopifnot(all.equal(get_cor(th0),th0/sqrt(1+th0^2)))
##' get_cor(c(0.5,0.2,0.5))
##' @export
get_cor <- function(theta) {
    n <- round((1  + sqrt(1+8*length(theta)))/2) ## dim of cor matrix
    L <- diag(n)
    L[lower.tri(L)] <- theta
    cL <- tcrossprod(L)
    Dh <- diag(1/sqrt(diag(cL)))
    cc <- Dh %*% cL %*% Dh
    return(cc[lower.tri(cc)])
}

match_which <- function(x,y) {
    which(sapply(y,function(z) x %in% z))
}

## reassign predvars to have term vars in the right order,
##  but with 'predvars' values inserted where appropriate
fix_predvars <- function(pv,tt) {
    if (length(tt)==3) {
        ## convert two-sided to one-sided formula
        tt <- RHSForm(tt, as.form=TRUE)
    }
    ## ugh, deparsing again ...
    tt_vars <- vapply(attr(tt, "variables"), deparse1, character(1))[-1]
    ## remove terminal paren - e.g. match term poly(x, 2) to
    ##   predvar poly(x, 2, <stuff>)
    ## beginning of string, including open-paren, colon
    ##  but not *first* comma nor arg ...
    ##  could possibly try init_regexp <- "^([^,]+).*" ?
    init_regexp <- "^([(^:_.[:alnum:]]+).*"
    tt_vars_short <- gsub(init_regexp,"\\1",tt_vars)
    if (is.null(pv) || length(tt_vars)==0) return(NULL)
    new_pv <- quote(list())
    ## maybe multiple variables per pv term ... [-1] ignores head
    ## FIXME: test for really long predvar strings ????
    pv_strings <- vapply(pv,deparse1,FUN.VALUE=character(1))[-1]
    pv_strings <- gsub(init_regexp,"\\1",pv_strings)
    for (i in seq_along(tt_vars)) {
        w <- match(tt_vars_short[[i]],pv_strings)
        if (!is.na(w)) {
            new_pv[[i+1]] <- pv[[w+1]]
        } else {
            ## insert symbol from term vars
            new_pv[[i+1]] <- as.symbol(tt_vars[[i]])
        }
    }
    return(new_pv)
}

hasRandom <- function(x) {
    pl <- getParList(x)
    return(length(unlist(pl[grep("^theta",names(pl))]))>0)
}

## retrieve parameters by name or index
getParms <- function(parm=NULL, object, full=FALSE) {
    vv <- vcov(object, full=TRUE)
    sds <- sqrt(diag(vv))
    pnames <- names(sds) <- rownames(vv)       ## parameter names (user-facing)
    intnames <- names(object$obj$env$last.par) ## internal names
    ## don't use object$obj$env$random; we want to keep "beta" vals, which may be
    ## counted as "random" if using REML
    intnames <- intnames[!intnames %in% c("b","bzi")]
    if (length(pnames) != length(sds)) { ## shouldn't happen ...
        stop("length mismatch between internal and external parameter names")
    }

    if (is.null(parm)) {
        if (!full && trivialDisp(object)) {
            parm <- grep("betad", intnames, invert=TRUE)
        } else {
            parm <- seq_along(sds)
        }
    }
    if (is.character(parm)) {
        if (identical(parm,"theta_")) {
            parm <- grep("^theta",intnames)
        } else if (identical(parm,"beta_")) {
            if (trivialDisp(object)) {
                ## include conditional and zi params
                ##   but not dispersion params
                parm <- grep("^beta(zi)?$",intnames)
            } else {
                parm <- grep("beta",intnames)
            }
        } else if (identical(parm, "disp_") ||
                   identical(parm, "sigma")) {
            parm <- grep("^betad", intnames)
        } else { ## generic parameter vector
            nparm <- match(parm,pnames)
            if (any(is.na(nparm))) {
                stop("unrecognized parameter names: ",
                     parm[is.na(nparm)])
            }
            parm <- nparm
        }
    }
    return(parm)
}

isREML <- function(x) {
    if (is.null(REML <- x$modelInfo$REML)) {
        ## let vcov work with old (pre-REML option) stored objects
        REML <- FALSE
    }
    return(REML)
}

## action: message, warning, stop
check_dots <- function(..., .ignore = NULL, .action="stop") {
    L <- list(...)
    if (length(.ignore)>0) {
        L <- L[!names(L) %in% .ignore]
    }
    if (length(L)>0) {
        FUN <- get(.action)
        FUN("unknown arguments: ",
            paste(names(L), collapse=","))
    }
    return(NULL)
}

if (getRversion()<"4.0.0") {
    deparse1 <- function (expr, collapse = " ", width.cutoff = 500L, ...) {
        paste(deparse(expr, width.cutoff, ...), collapse = collapse)
    }
}

## in case these are useful, we can document and export them later ...
#' @importFrom stats rnbinom qnbinom dnbinom pnbinom

rnbinom1 <- function(n, mu, phi) {
    ## var = mu*(1+phi) = mu*(1+mu/k) -> k = mu/phi
    rnbinom(n, mu=mu, size=mu/phi)
}

dnbinom1 <- function(x, mu, phi, ...) {
    dnbinom(x, mu=mu, size=mu/phi, ...)
}

pnbinom1 <- function(q, mu, phi, ...) {
    pnbinom(q, mu=mu, size=mu/phi, ...)
}

qnbinom1 <- function(p, mu, phi, ...) {
    qnbinom(p, mu=mu, size=mu/phi, ...)
}

nullSparseMatrix <- function() {
    argList <- list(
        dims=c(0,0),
        i=integer(0),
        j=integer(0),
        x=numeric(0))
    if (utils::packageVersion("Matrix")<"1.3.0") {
        do.call(Matrix::sparseMatrix, c(argList, list(giveCsparse=FALSE)))
    } else {
        do.call(Matrix::sparseMatrix, c(argList, list(repr="T")))
    }
}


##' Check OpenMP status
##'
##' Checks whether OpenMP has been successfully enabled for this
##' installation of the package. (Use the \code{parallel} argument
##' to \code{\link{glmmTMBControl}}, or set \code{options(glmmTMB.cores=[value])},
##' to specify that computations should be done in parallel.)
##' @seealso \code{\link[TMB]{benchmark}}, \code{\link{glmmTMBControl}}
##' @return \code{TRUE} or {FALSE} depending on availability of OpenMP
##' @export
omp_check <- function() {
    .Call("omp_check", PACKAGE="glmmTMB")
}

get_pars <- function(object, unlist=TRUE) {
    ee <- object$obj$env
    x <- ee$last.par.best
    ## work around built-in default to parList, which
    ##  is bad if no random component
    if (length(ee$random)>0) x <- x[-ee$random]
    p <- ee$parList(x=x)
    if (!unlist) return(p)
    p <- unlist(p[names(p)!="b"])  ## drop primary RE
    names(p) <- gsub("[0-9]+$","",names(p)) ## remove disambiguators
    return(p)
}

