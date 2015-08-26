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

##' @importFrom stats reformulate
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

##' combine unary or binary operator + arguments (sugar for 'substitute')
makeOp <- function(x,y,op=NULL) {
    if (is.null(op)) {  ## unary
        substitute(OP(X),list(X=x,OP=y))
    } else substitute(OP(X,Y), list(X=x,OP=op,Y=y))
}

##' combines the right-hand sides of two formulas
##' @param f1
##' @param f2
##' @keywords internal
##' @examples
##' if(FALSE) ## not exported
##' addForm0(y~x,~1)
##' addForm0(~x,~y)
addForm0 <- function(f1,f2) {
  if (length(f2)==3) warning("discarding RHS of second argument")
  RHSForm(f1) <- makeOp(RHSForm(f1),RHSForm(f2),quote(`+`))
  return(f1)
}

##' combine right-hand sides of an arbitrary number of formulas
addForm <- function(...) {
  Reduce(addForm0,list(...))
}

##' deparse(.) returning \bold{one} string
##' @note Protects against the possibility that results from deparse() will be
##'       split after 'width.cutoff' (by default 60, maximally 500)
safeDeparse <- function(x, collapse=" ") {
    paste(deparse(x, 500L), collapse=collapse)
}

##' list of specials -- taken from enum.R
findReTrmClasses <- function() {
    names(.valid_covstruct)
}

##' find individual slashless terms within
##' formula elements
slashTerms <- function(x,debug=FALSE)
{
    if (debug) cat("slashTerms: ",deparse(x),"\n")
    if (!("/" %in% all.names(x))) return(x)
    if (x[[1]] == as.name("("))  ## strip parens if necessary
        return(slashTerms(x[[2]]))
    if (x[[1]] != as.name("/")) {
        stop("didn't expect head ",deparse(x[[1]]))
    }
    ## recursively apply to all terms
    lapply(x[-1],slashTerms,debug=debug)
}


##' Expand any slashes in the grouping factors returned by fb
##' 
##' Original by Doug Bates: copied and streamlined from lme4/utilities.R
##' FIXME: expand to allow expansion of '+', '*'
##' @examples
##' ff <- fbx(y~1+(x|f/g))
##' expandSlash(ff)
##' expandSlash(quote(1|(f/g)/h))
##' expandSlash(quote(1|f/g/h))
expandSlash <- function(bb,debug=FALSE) {
    ## Create the interaction terms for nested effects
    makeInteraction <- function(x) {
        if (length(x) < 2) return(x)
        trm1 <- makeInteraction(x[[1]])
        trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
        list(makeOp(x[[2]],trm11,quote(`:`)), trm1)
    }
    ## Return the list of '/'-separated terms
    if (!is.list(bb))
        expandSlash(list(bb),debug=debug)
    else {
        for (i in seq_along(bb)) {
            esfun <- function(x) {
                if (!("/" %in% all.names(x))) return(x)
                if (length(x)==1) stop("term should have length >1")
                if (length(x)==2) {
                    ## unary operator such as diag(1|f/g)
                    ## return diag(...) + diag(...) + ...
                    return(lapply(esfun(x[[2]]),
                                  makeOp,y=head(x)))
                }
                if (!head(x)==quote(`|`)) stop("operator should be a bar")
                if (is.list(trms <- slashTerms(x[[3]]))) {
                    ## replicate LHS of bar with each RHS term
                    return(lapply(unlist(makeInteraction(trms)),
                                  makeOp,x=x[[2]],op=quote(`|`)))
                }  else {
                    return(x)
                }
            } ## esfun def.
            return(unlist(lapply(bb,esfun)))
        } ## loop over bb
    }
} ## {expandSlash}

## sugar
head.formula <- head.call <- function(x, ...) {
    x[[1]]
}

##' (f)ind (b)ars e(x)tended: recursive
##' 
##' 1. atom (not a call or an expression): NULL
##' 2. special, i.e. foo(...) where "foo" is in specials: return term
##' 3. parenthesized term: *if* the head of the head is | (i.e.
##'    it is of the form (xx|gg), then convert it to the default
##'    special type; we won't allow pathological cases like
##'    ((xx|gg)) ... [can we detect them?]
fbx <- function(term,debug=FALSE,specials=character(0),
                default.special="us") {
    ds <- eval(substitute(as.name(foo),list(foo=default.special)))
    if (is.name(term) || !is.language(term)) return(NULL)
    if (list(term[[1]]) %in% lapply(specials,as.name)) {
        if (debug) cat("special: ",deparse(term),"\n")
        return(term)
    }
    if (head(term) == as.name('|')) {  ## found x | g
        if (debug) cat("bar term:",deparse(term),"\n")
        return(makeOp(term,ds))
    }
    if (head(term) == as.name("(")) {  ## found (...)
        if (debug) cat("paren term:",deparse(term),"\n")
        return(fbx(term[[2]],debug,specials))
    }
    stopifnot(is.call(term))
    if (length(term) == 2) {
        ## unary operator, decompose argument
        if (debug) cat("unary operator:",deparse(term[[2]]),"\n")
        return(fbx(term[[2]],debug,specials))
    }
    ## binary operator, decompose both arguments
    if (debug) cat("binary operator:",deparse(term[[2]]),",",
                   deparse(term[[3]]),"\n")
    c(fbx(term[[2]],debug,specials), fbx(term[[3]],debug,specials))
}

##' Parse a formula into fixed formula and random effect terms,
##' treating 'special' terms (of the form foo(x|g[,m])) appropriately
##'
##' Taken from Steve Walker's lme4ord,
##' ultimately from the flexLambda branch of lme4
##' <https://github.com/stevencarlislewalker/lme4ord/blob/master/R/formulaParsing.R>.  Mostly for internal use.
##' @title Split formula containing special random effect terms
##' @param formula a formula containing special random effect terms
##' @param defaultTerm default type for non-special RE terms
##' @param allowFixedOnly (logical) are formulas with no RE terms OK?
##' @param allowNoSpecials (logical) are formulas with only standard RE terms OK?
##' @return a list containing elements \code{fixedFormula};
##' \code{reTrmFormulas} list of \code{x | g} formulas for each term;
##' \code{reTrmAddArgs} list of function+additional arguments, i.e. \code{list()} (non-special), \code{foo()} (no additional arguments), \code{foo(addArgs)} (additional arguments); \code{reTrmClasses} (vector of special functions/classes, as character)
##' @examples
##' splitForm(~x+y)                     ## no specials or RE
##' splitForm(~x+y+(f|g))               ## no specials
##' splitForm(~x+y+diag(f|g))           ## one special
##' splitForm(~x+y+(diag(f|g)))         ## 'hidden' special
##' splitForm(~x+y+(f|g)+cs(1|g))       ## combination
##' splitForm(~x+y+(1|f/g))             ## 'slash'; term
##' splitForm(~x+y+(1|f/g/h))             ## 'slash'; term
##' splitForm(~x+y+(1|(f/g)/h))             ## 'slash'; term
##' splitForm(~x+y+(f|g)+cs(1|g)+cs(a|b,stuff))  ## complex special
##' splitForm(~(((x+y))))               ## lots of parentheses
##' 
##' @author Steve Walker
##' @importFrom lme4 nobars
##' @export
splitForm <- function(formula,
                      defaultTerm="us",
                      allowFixedOnly=TRUE,
                      allowNoSpecials=TRUE,
                      debug=FALSE) {

    ## logic:
    
    ## string for error message *if* specials not allowed
    ## (probably package-specific)
    noSpecialsAlt <- "lmer or glmer"

    specials <- findReTrmClasses()
    
    ## formula <- expandDoubleVerts(formula)
    ## split formula into separate
    ## random effects terms
    ## (including special terms)

    formSplits <- fbx(formula,debug,specials)
    formSplits <- expandSlash(formSplits,debug)
                                        # check for hidden specials
                                        # (i.e. specials hidden behind
                                        # parentheses)
    if (length(formSplits)>0) {
        formSplitID <- sapply(lapply(formSplits, "[[", 1), as.character)
                                        # warn about terms without a
                                        # setReTrm method

        ## FIXME:: do we need all of this??

        if (FALSE) {
            badTrms <- formSplitID == "|"
        ## if(any(badTrms)) {
        ## stop("can't find setReTrm method(s)\n",
        ## "use findReTrmClasses() for available methods")
        ## FIXME: coerce bad terms to default as attempted below
        ## warning(paste("can't find setReTrm method(s) for term number(s)",
        ## paste(which(badTrms), collapse = ", "),
        ## "\ntreating those terms as unstructured"))
            formSplitID[badTrms] <- "("
            fixBadTrm <- function(formSplit) {
                makeOp(formSplit[[1]],quote(`(`))
                ## as.formula(paste(c("~(", as.character(formSplit)[c(2, 1, 3)], ")"),
                ## collapse = " "))[[2]]
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

    list(fixedFormula  = fixedFormula,
         reTrmFormulas = reTrmFormulas,
         reTrmAddArgs  = reTrmAddArgs,
         reTrmClasses  = reTrmClasses)
}

##' @param term language object
##' @rdname splitForm
##' @examples
##' noSpecials(y~1+us(1|f))
##' noSpecials(y~1+us(1|f),delete=FALSE)
##' @export
##' @keywords internal
noSpecials <- function(term, delete=TRUE) {
    nospec <- noSpecials_(term, delete=delete)
    if (inherits(term, "formula") && length(term) == 3 && is.symbol(nospec)) {
        ## called with two-sided RE-only formula:
        ##    construct response~1 formula
        reformulate("1", response = deparse(nospec))
    } else
        nospec
}


## noSpecials_(y~1+us(1|f))
noSpecials_ <- function(term,delete=TRUE) {
    if (!anySpecial(term)) return(term)
    if (isSpecial(term)) {
        if(delete) {
            return(NULL)
        } else {
            ## FIXME: returns 1 | f, would like (1|f)
            return(term[[2]])
        }
    }
    nb2 <- noSpecials(term[[2]],delete=delete)
    nb3 <- noSpecials(term[[3]],delete=delete)
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
    for(tt in term)
        if(isSpecial(tt)) return(TRUE)
    FALSE
}

## FIXME: this could be fooled by a term with a matching name
## should really look for [special]\\(.+\\)
anySpecial <- function(term) {
    any(findReTrmClasses() %in% all.names(term))
}
