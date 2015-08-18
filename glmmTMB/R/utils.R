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
