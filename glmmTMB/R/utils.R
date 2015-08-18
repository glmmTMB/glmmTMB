## generate a list with names equal to values
namedList <- function (...) 
{
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L))) 
        nm <- snm
    if (any(nonames <- nm == "")) 
        nm[nonames] <- snm[nonames]
    setNames(L, nm)
}

RHSForm <- function(form,as.form=FALSE) {
    rhsf <- form[[length(form)]]
    if (as.form) reformulate(deparse(rhsf)) else rhsf
}

`RHSForm<-` <- function(formula,value) {
    formula[[length(formula)]] <- value
    formula
}

##' Random Effects formula only
reOnly <- function(f,response=FALSE) {
    response <- if (response && length(f)==3) f[[2]] else NULL
    reformulate(paste0("(", vapply(findbars(f), safeDeparse, ""), ")"),
                response=response)
}

##' deparse(.) returning \bold{one} string
##' @note Protects against the possibility that results from deparse() will be
##'       split after 'width.cutoff' (by default 60, maximally 500)
safeDeparse <- function(x, collapse=" ") paste(deparse(x, 500L), collapse=collapse)


##' .. content for \description{} (no empty lines) ..
##'
##' Stolen from flexLambda branch 
##' @title split formula containing 'specials'
##' @param formula 
##' @param specials 
##' @return something
##' @examples
##' glmmTMB:::splitForm(~x+y+diag(f|g))
##' @author Fabian Scheipl?
##' @importFrom lme4 nobars
##' @keywords internal
splitForm <- function(formula, specials = c("diag", "cs", "ar1d")) {

    ## Recursive function: (f)ind (b)ars (a)nd (s)pecials
    ## cf. fb function in findbars (i.e. this is a little DRY)
    fbas <- function(term) {
        if (is.name(term) || !is.language(term)) return(NULL)
        for (sp in specials) if (term[[1]] == as.name(sp)) return(term)
        if (term[[1]] == as.name("(")) return(term)
        stopifnot(is.call(term))
        if (term[[1]] == as.name('|')) return(term)
        if (length(term) == 2) return(fbas(term[[2]]))
        c(fbas(term[[2]]), fbas(term[[3]]))
    }
    ## we don't need this
    ## formula <- expandDoubleVerts(formula)
                                        # split formula into separate
                                        # random effects terms
                                        # (including special terms)
    formSplits <- fbas(formula)
                                        # vector to identify what
                                        # special (by name), or give
                                        # "(" for standard terms
    formSplitID <- sapply(lapply(formSplits, "[[", 1), as.character)
                                        # standard RE terms
    formSplitStan <- formSplits[formSplitID == "("]
                                        # special RE terms
    formSplitSpec <- formSplits[!(formSplitID == "(")]

    if(length(formSplitSpec) == 0) stop(
             "no special covariance structures. ",
             "please use lmer, not flexLmer")

                                        # construct the formula for
                                        # the specials only
    reGenerators <- as.formula(paste("~ ", paste(formSplitSpec, collapse = " + ")))   


                                        # construct the standard
                                        # 'lmer'-style formula
                                        # (without the specials)
    if(length(formSplitStan) == 0) {
        lmerformula <- nobars(formula)
    }
    return(lmerformula)
}
