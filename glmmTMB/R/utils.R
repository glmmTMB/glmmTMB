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
safeDeparse <- function(x, collapse=" ") {
    paste(deparse(x, 500L), collapse=collapse)
}

##' list of specials
##' FIXME: 
findReTrmClasses <- function() {
    c("diag","cs","us")
}

##' Parse a formula into fixed formula and random effect terms,
##' treating 'special' terms (of the form foo(x|g[,m])) appropriately
##'
##' Taken from Steve Walker's lme4ord,
##' ultimately from Fabian Scheipl's flexLambda branch of lme4
##' <https://github.com/stevencarlislewalker/lme4ord/blob/master/R/formulaParsing.R>
##' @title Split formula containing special random effect terms
##' @param formula a formula containing special random effect terms
##' @param defaultTerm default type for non-special RE terms
##' @param addSpecials additional special types not in standard list
##' @param allowFixedOnly (logical) are formulas with no RE terms OK?
##' @param allowNoSpecials (logical) are formulas with only standard RE terms OK?
##' @param noSpecialsAlt name for alternative functions to be used with non-special formulas
##' @return a list containing elements
##' @examples
##' glmmTMB:::splitForm(~x+y)            ## no specials or RE
##' glmmTMB:::splitForm(~x+y+(f|g))      ## no specials
##' glmmTMB:::splitForm(~x+y+diag(f|g))  ## one special
##' glmmTMB:::splitForm(~x+y+(f|g)+cs(1|g))
##' glmmTMB:::splitForm(~x+y+(f|g)+cs(1|g)+foo(a|b,stuff),
##'                     addSpecials=c("cs","foo"))
##' @author Fabian Scheipl, Steve Walker
##' @importFrom lme4 nobars
##' @export 
splitForm <- function(formula,defaultTerm="unstruc",
                      addSpecials=NULL,
                      allowFixedOnly=TRUE,
                      allowNoSpecials=TRUE,
                      noSpecialsAlt="lmer or glmer") {
    
    specials <- c(findReTrmClasses(),addSpecials)
    ## ignore any specials not in formula
    specialsToKeep <- vapply(specials, grepl,
                             x = safeDeparse(formula),
                             logical(1))
    specials <- specials[specialsToKeep]

    ## Recursive function: (f)ind (b)ars (a)nd (s)pecials
    ## cf. fb function in findbars (i.e. this is a little DRY)
    fbas <- function(term,debug=FALSE) {
        if (is.name(term) || !is.language(term)) return(NULL)
        for (sp in specials) if (term[[1]] == as.name(sp)) { ## found special
            if (debug) cat("special: ",deparse(term),"\n")
            return(term)
        }
        if (term[[1]] == as.name("(")) {  ## found (...)
            if (debug) cat("paren term:",deparse(term),"\n")
            return(term)
        }
        stopifnot(is.call(term))
        if (term[[1]] == as.name('|')) {  ## found x | g
            if (debug) cat("bar term:",deparse(term),"\n")
            return(term)
        }
        if (length(term) == 2) {
            ## unary operator, decompose argument
            if (debug) cat("unary operator:",deparse(term[[2]]),"\n")
            return(fbas(term[[2]],debug=debug))
        }
        ## binary operator, decompose both arguments
        if (debug) cat("binary operator:",deparse(term[[2]]),",",
                       deparse(term[[3]]),"\n")
        c(fbas(term[[2]],debug=debug), fbas(term[[3]],debug=debug))
    }

    ## not used in glmmTMB
    ## formula <- expandDoubleVerts(formula)
                                        # split formula into separate
                                        # random effects terms
                                        # (including special terms)
    formSplits <- fbas(formula,debug=TRUE)
                                        # check for hidden specials
                                        # (i.e. specials hidden behind
                                        # parentheses)

    if (!is.null(formSplits)) {
        formSplits <- lapply(formSplits, uncoverHiddenSpecials)
                                        # vector to identify what
                                        # special (by name), or give
                                        # "(" for standard terms, or
                                        # give "|" for specials
                                        # without a setReTrm method
        formSplitID <- sapply(lapply(formSplits, "[[", 1), as.character)
        as.character(formSplits[[1]])
                                        # warn about terms without a
                                        # setReTrm method
        badTrms <- formSplitID == "|"
        if(any(badTrms)) {
            stop("can't find setReTrm method(s)\n",
                 "use findReTrmClasses() for available methods")
            ## FIXME: coerce bad terms to default as attempted below
            warning(paste("can't find setReTrm method(s) for term number(s)",
                          paste(which(badTrms), collapse = ", "),
                          "\ntreating those terms as unstructured"))
            formSplitID[badTrms] <- "("
            fixBadTrm <- function(formSplit) {
                as.formula(paste(c("~(", as.character(formSplit)[c(2, 1, 3)], ")"),
                                 collapse = " "))[[2]]
            }
            formSplits[badTrms] <- lapply(formSplits[badTrms], fixBadTrm)
        }

        parenTerm <- formSplitID == "("
                                        # capture additional arguments
        reTrmAddArgs <- lapply(formSplits, "[", -2)[!parenTerm]
                                        # remove these additional
                                        # arguments
        formSplits <- lapply(formSplits, "[", 1:2)
                                        # standard RE terms
        formSplitStan <- formSplits[parenTerm]
                                        # structured RE terms
        formSplitSpec <- formSplits[!parenTerm]

        if (!allowNoSpecials) {
            if(length(formSplitSpec) == 0) stop(
                     "no special covariance structures. ",
                     "please use ",noSpecialsAlt,
                     " or use findReTrmClasses() for available structures.")
        }
        
        reTrmFormulas <- c(lapply(formSplitStan, "[[", 2),
                           lapply(formSplitSpec, "[[", 2))
        reTrmClasses <- c(rep(defaultTerm, length(formSplitStan)),
                          sapply(lapply(formSplitSpec, "[[", 1), as.character))
    } else {
        reTrmFormulas <- reTrmAddArgs <- reTrmClasses <- NULL
    }
    fixedFormula <- noSpecials(nobars(formula))


    return(list(fixedFormula  = fixedFormula,
                reTrmFormulas = reTrmFormulas,
                reTrmAddArgs  = reTrmAddArgs,
                reTrmClasses  = reTrmClasses))
}

uncoverHiddenSpecials <- function(trm) {
    if(trm[[1]] == "(") {
        if(anySpecial(trm[[2]][[1]])) trm <- trm[[2]]
    }
    return(trm)
}

##' @param term language object
##' @rdname splitForm
##' @keywords internal
noSpecials <- function(term) {
    nospec <- noSpecials_(term)
    if (is(term,"formula") && length(term)==3 && is.symbol(nospec)) {
        ## called with two-sided RE-only formula:
        ##    construct response~1 formula
        nospec <- reformulate("1", response = deparse(nospec))
    }
    return(nospec)
}

noSpecials_ <- function(term) {
    if (!anySpecial(term)) return(term)
    if (isSpecial(term)) return(NULL)
    nb2 <- noSpecials(term[[2]])
    nb3 <- noSpecials(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

isSpecial <- function(term) {
    if(is.call(term)) {
        for(cls in findReTrmClasses()) {
            if(term[[1]] == cls) return(TRUE)
        }
    }
    FALSE
}

isAnyArgSpecial <- function(term) {
    for(i in seq_along(term)) {
        if(isSpecial(term[[i]])) return(TRUE)
    }
    FALSE
}

anySpecial <- function(term) {
    any(findReTrmClasses() %in% all.names(term))
}
