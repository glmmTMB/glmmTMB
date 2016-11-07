## backward compat (copied from lme4)
if((Rv <- getRversion()) < "3.2.1") {
    lengths <- function (x, use.names = TRUE) vapply(x, length, 1L, USE.NAMES = use.names)
}
rm(Rv)

## generate a list with names equal to values
namedList <- function (...) {
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
    if (as.form) as.formula(substitute(~F,list(F=rhsf)),
                            env=environment(form)) else rhsf
}

`RHSForm<-` <- function(formula,value) {
    formula[[length(formula)]] <- value
    formula
}

## Random Effects formula only
## reOnly <- function(f,response=FALSE) {
##    response <- if (response && length(f)==3) f[[2]] else NULL
##    reformulate(paste0("(", vapply(findbars(f), safeDeparse, ""), ")"),
##                response=response)
## }

## better version -- operates on language objects (no deparse())
reOnly <- function(f,response=FALSE,bracket=TRUE) {
    ff <- f
    if (bracket)
        ff <- lapply(findbars(ff),makeOp,quote(`(`)) ## bracket-protect terms
    ff <- Reduce(function(x,y) makeOp(x,y,op=quote(`+`)),ff)
    if (response && length(f)==3) {
        form <- makeOp(f[[2]],ff,quote(`~`))
    } else {
        form <- makeOp(ff,quote(`~`))
    }
    return(form)
}

## combine unary or binary operator + arguments (sugar for 'substitute')
## FIXME: would be nice to have multiple dispatch, so
## (arg,op) gave unary, (arg,arg,op) gave binary operator
makeOp <- function(x,y,op=NULL) {
    if (is.null(op)) {  ## unary
        substitute(OP(X),list(X=x,OP=y))
    } else substitute(OP(X,Y), list(X=x,OP=op,Y=y))
}

## combines the right-hand sides of two formulas
## @param f1 formula #1
## @param f2 formula #2
## @examples
## if (FALSE) {  ## still being exported despite "<at>keywords internal" ??
## addForm0(y~x,~1)
## addForm0(~x,~y)
## }
## @keywords internal
addForm0 <- function(f1,f2) {
  if (length(f2)==3) warning("discarding LHS of second argument")
  RHSForm(f1) <- makeOp(RHSForm(f1),RHSForm(f2),quote(`+`))
  return(f1)
}

##' Combine right-hand sides of an arbitrary number of formulas
##' @param ... arguments to pass through to \code{addForm0}
##' @rdname splitForm
##' @export
addForm <- function(...) {
  Reduce(addForm0,list(...))
}

addArgs <- function(argList) {
  Reduce(function(x,y) makeOp(x,y,op=quote(`+`)),argList)
}

##' list of specials -- taken from enum.R
findReTrmClasses <- function() {
    names(.valid_covstruct)
}

## expandGrpVar(quote(x*y))
## expandGrpVar(quote(x/y))
expandGrpVar <- function(f) {
    form <- as.formula(makeOp(f,quote(`~`)))
    mm <- terms(form)
    toLang <- function(x) parse(text=x)[[1]]
    lapply(attr(mm,"term.labels"),
           toLang)
}

##' expand interactions/combinations of grouping variables
##'
##' Modeled after lme4:::expandSlash, by Doug Bates
##' @param bb a list of naked grouping variables, i.e. 1 | f
##' @examples
##' ff <- glmmTMB:::fbx(y~1+(x|f/g))
##' glmmTMB:::expandAllGrpVar(ff)
##' glmmTMB:::expandAllGrpVar(quote(1|(f/g)/h))
##' glmmTMB:::expandAllGrpVar(quote(1|f/g/h))
##' glmmTMB:::expandAllGrpVar(quote(1|f*g))
##' @importFrom utils head
##' @keywords internal
expandAllGrpVar <- function(bb) {
        ## Return the list of '/'-separated terms
    if (!is.list(bb))
        expandAllGrpVar(list(bb))
    else {
        for (i in seq_along(bb)) {
            esfun <- function(x) {
                if (length(x)==1) return(x)
                if (length(x)==2) {
                    ## unary operator such as diag(1|f/g)
                    ## return diag(...) + diag(...) + ...
                    return(lapply(esfun(x[[2]]),
                                  makeOp,y=head(x)))
                }
                if (length(x)==3) {
                    ## binary operator
                    if (x[[1]]==quote(`|`)) {
                        return(lapply(expandGrpVar(x[[3]]),
                                      makeOp,x=x[[2]],op=quote(`|`)))
                    } else {
                        return(makeOp(esfun(x[[2]]),esfun(x[[3]]),
                                      op=x[[1]]))
                    }
                }
            } ## esfun def.
            return(unlist(lapply(bb,esfun)))
        } ## loop over bb
    }
}

## sugar -- ???? this returns '~'  ???
head.formula <- head.call <- function(x, ...) {
    x[[1]]
}

##' (f)ind (b)ars e(x)tended: recursive
##'
##' @param term a formula or piece of a formula
##' @param debug (logical) debug?
##' @param specials list of special terms
##' @param default.special character: special to use for parenthesized terms - i.e. random effects terms with unspecified structure
##' 1. atom (not a call or an expression): NULL
##' 2. special, i.e. foo(...) where "foo" is in specials: return term
##' 3. parenthesized term: \emph{if} the head of the head is | (i.e.
##'    it is of the form (xx|gg), then convert it to the default
##'    special type; we won't allow pathological cases like
##'    ((xx|gg)) ... [can we detect them?]
##' @keywords internal
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

    fbxx <- fbx(formula,debug,specials)
    formSplits <- expandAllGrpVar(fbxx)

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

        }  ## skipped

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
##' noSpecials(y~us(1|f))
##' @export
##' @keywords internal
noSpecials <- function(term, delete=TRUE) {
    nospec <- noSpecials_(term, delete=delete)
    if (inherits(term, "formula") && length(term) == 3 && is.symbol(nospec)) {
        ## called with two-sided RE-only formula:
        ##    construct response~1 formula
        as.formula(substitute(R~1,list(R=nospec)),
                   env=environment(term))
    } else
        nospec
}


## noSpecials_(y~1+us(1|f))
noSpecials_ <- function(term,delete=TRUE) {
    if (!anySpecial(term)) return(term)
    if (isSpecial(term)) {
        if(delete) {
            NULL
        } else { ## careful to return  (1|f) and not  1|f:
            substitute((TERM), list(TERM = term[[2]]))
        }
    } else {
        nb2 <- noSpecials(term[[2]],delete=delete)
        nb3 <- if (length(term)==3) {
                   noSpecials(term[[3]],delete=delete)
               } else NULL
        if (is.null(nb2))
            nb3
        else if (is.null(nb3))
            nb2
        else {
            term[[2]] <- nb2
            term[[3]] <- nb3
            term
        }
    }
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

## does the formula contain a particular value?
## inForm(z~.,quote(.))
## inForm(z~y,quote(.))
## inForm(z~a+b+c,quote(c))
## inForm(z~a+b+(d+e),quote(c))
inForm <- function(form,value) {
    if (any(sapply(form,identical,value))) return(TRUE)
    if (all(sapply(form,length)==1)) return(FALSE)
    return(any(sapply(form,inForm,value)))
}

# drop.special(x~a + b+ offset(z))
drop.special <- function(term,value=quote(offset)) {
    if (length(term)==2 && identical(term[[1]],value)) return(NULL)
    if (length(term)==1) return(term)
    ## recurse, treating unary and binary operators separately
    nb2 <- drop.special(term[[2]])
    nb3 <- if (length(term)==3) {
               drop.special(term[[3]])
           } else NULL
    if (is.null(nb2)) ## RHS was special-only
        nb3
    else if (is.null(nb3)) ## LHS was special-only
        nb2
    else {
        ## insert values into daughters and return
        term[[2]] <- nb2
        term[[3]] <- nb3
        return(term)
    }
}

## from Gabor Grothendieck: recursive solution
## http://stackoverflow.com/questions/40308944/removing-offset-terms-from-a-formula
drop.special2 <- function(x, value=quote(offset), preserve = NULL) {
  k <- 0
  proc <- function(x) {
    if (length(x) == 1) return(x)
    if (x[[1]] == value && !((k<<-k+1) %in% preserve)) return(x[[1]])
    replace(x, -1, lapply(x[-1], proc))
  }
  update(proc(x), substitute(. ~ . - x,list(x=value)))
}
