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

sumTerms <- function(termList) {
    Reduce(function(x,y) makeOp(x,y,op=quote(`+`)),termList)
}

## better version -- operates on language objects (no deparse())
reOnly <- function(f,response=FALSE,bracket=TRUE) {
    ff <- f
    if (bracket)
        ff <- lapply(findbars(ff),makeOp,quote(`(`)) ## bracket-protect terms
    ff <- sumTerms(ff)
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
    if (is.null(op) || missing(y)) {  ## unary
        if (is.null(op)) {
            substitute(OP(X),list(X=x,OP=y))
        } else {
            substitute(OP(X),list(X=x,OP=op))
        }
    } else substitute(OP(X,Y), list(X=x,OP=op,Y=y))
}

## combines the right-hand sides of two formulas, or a formula and a symbol
## @param f1 formula #1
## @param f2 formula #2
## @examples
## if (FALSE) {  ## still being exported despite "<at>keywords internal" ??
## addForm0(y~x,~1)
## addForm0(~x,~y)
## }
## @keywords internal
addForm0 <- function(f1,f2,naked=FALSE) {
    tilde <- as.symbol("~")
    if (!identical(head(f2),tilde)) {
        f2 <- makeOp(f2,tilde)
    }
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

## sugar: this returns the operator, whether ~ or something else
head.formula <- head.call <- function(x, ...) {
    x[[1]]
}

## sugar: we can call head on a symbol and get back the symbol
head.name <- function(x) { x }

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
##' @param debug debugging mode (print stuff)?
##' @examples
##' noSpecials(y~1+us(1|f))
##' noSpecials(y~1+us(1|f),delete=FALSE)
##' noSpecials(y~us(1|f))
##' noSpecials(y~us+1)  ## should *not* delete unless head of a function
##' @export
##' @keywords internal
noSpecials <- function(term, delete=TRUE, debug=FALSE) {
    nospec <- noSpecials_(term, delete=delete, debug=debug)
    if (inherits(term, "formula") && length(term) == 3 && is.symbol(nospec)) {
        ## called with two-sided RE-only formula:
        ##    construct response~1 formula
        as.formula(substitute(R~1,list(R=nospec)),
                   env=environment(term))
    } else
        nospec
}


## noSpecials_(y~1+us(1|f))
noSpecials_ <- function(term,delete=TRUE, debug=FALSE) {
    if (debug) print(term)
    if (!anySpecial(term)) return(term)
    if (length(term)==1) return(term)  ## 'naked' specials
    if (isSpecial(term)) {
        if(delete) {
            NULL
        } else { ## careful to return  (1|f) and not  1|f:
            substitute((TERM), list(TERM = term[[2]]))
        }
    } else {
        nb2 <- noSpecials(term[[2]], delete=delete, debug=debug)
        nb3 <- if (length(term)==3) {
                   noSpecials(term[[3]], delete=delete, debug=debug)
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
        ## %in% doesn't work (requires vector args)
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

## This could be in principle be fooled by a term with a matching name
## but this case is caught in noSpecials_() where we test for length>1
anySpecial <- function(term) {
    any(findReTrmClasses() %in% all.names(term))
}

##' test formula: does it contain a particular element?
##' @rdname formFuns
##' @examples
##' inForm(z~.,quote(.))
##' inForm(z~y,quote(.))
##' inForm(z~a+b+c,quote(c))
##' inForm(z~a+b+(d+e),quote(c))
##' f <- ~ a + offset(x)
##' f2 <- z ~ a
##' inForm(f,quote(offset))
##' inForm(f2,quote(offset))
##' @export
##' @keywords internal
inForm <- function(form,value) {
    if (any(sapply(form,identical,value))) return(TRUE)
    if (all(sapply(form,length)==1)) return(FALSE)
    return(any(vapply(form,inForm,value,FUN.VALUE=logical(1))))
}

##' extract terms with a given head from an expression/formula
##' @rdname formFuns
##' @param term expression/formula
##' @param value head of terms to extract
##' @return a list of expressions
##' @examples
##' extractForm(~a+offset(b),quote(offset))
##' extractForm(~c,quote(offset))
##' extractForm(~a+offset(b)+offset(c),quote(offset))
##' @export
##' @keywords internal
extractForm <- function(term,value) {
    if (!inForm(term,value)) return(NULL)
    if (is.name(term) || !is.language(term)) return(NULL)
    if (identical(head(term),value)) {
        return(term)
    }
    if (length(term) == 2) {
        return(extractForm(term[[2]],value))
    }
    return(c(extractForm(term[[2]],value),
             extractForm(term[[3]],value)))
}

##' return a formula/expression with a given value stripped, where
##' it occurs as the head of a term
##' @rdname formFuns
##' @examples 
##' dropHead(~a+offset(b),quote(offset))
##' dropHead(~a+poly(x+z,3)+offset(b),quote(offset))
##' @export
##' @keywords internal
dropHead <- function(term,value) {
    if (!inForm(term,value)) return(term)
    if (is.name(term) || !is.language(term)) return(term)
    if (identical(head(term),value)) {
        return(term[[2]])
    }
    if (length(term) == 2) {
        return(dropHead(term[[2]],value))
    } else  if (length(term) == 3) {
        term[[2]] <- dropHead(term[[2]],value)
        term[[3]] <- dropHead(term[[3]],value)
        return(term)
    } else stop("length(term)>3")
}


## UNUSED (same function as drop.special2?)
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

##' drop terms matching a particular value from an expression
##' @rdname formFuns
## from Gabor Grothendieck: recursive solution
## http://stackoverflow.com/questions/40308944/removing-offset-terms-from-a-formula
##' @param x formula
##' @param value term to remove from formula
##' @param preserve (integer) retain the specified occurrence of "value"
##' @keywords internal
drop.special2 <- function(x, value=quote(offset), preserve = NULL) {
  k <- 0
  proc <- function(x) {
    if (length(x) == 1) return(x)
    if (x[[1]] == value && !((k <<- k+1) %in% preserve)) return(x[[1]])
    replace(x, -1, lapply(x[-1], proc))
  }
  ## handle 1- and 2-sided formulas
  if (length(x)==2) {
      newform <- substitute(~ . -x, list(x=value))
  } else {
      newform <- substitute(. ~ . - x, list(x=value))
  }
  return(update(proc(x), newform))
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
