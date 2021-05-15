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

RHSForm <- function(form,as.form=FALSE) {
    if (!as.form) return(form[[length(form)]])
    if (length(form)==2) return(form)  ## already RHS-only
    ## by operating on RHS in situ rather than making a new formula
    ## object, we avoid messing up existing attributes/environments etc.
    form[[2]] <- NULL
    ## assumes response is *first* variable (I think this is safe ...)
    if (length(vars <- attr(form,"variables"))>0) {
        attr(form,"variables") <- vars[-2]
    }
    if (is.null(attr(form,"response"))) {
        attr(form,"response") <- 0
    }
    if (length(facs <- attr(form,"factors"))>0) {
        attr(form,"factors") <- facs[-1,]
    }
    return(form)
}

`RHSForm<-` <- function(formula,value) {
    formula[[length(formula)]] <- value
    formula
}


sumTerms <- function(termList) {
    Reduce(function(x,y) makeOp(x,y,op=quote(`+`)),termList)
}

## extract random effects component of formula
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
                        return(setNames(makeOp(esfun(x[[2]]),esfun(x[[3]]),
                                               op=x[[1]]),names(x)))
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
##' @examples
##' splitForm(quote(us(x,n=2)))
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
##' splitForm(~1+rr(f|g,n=2))
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
        reTrmFormulas <- unlist(reTrmFormulas) # Fix me:: added for rr structure when it has n = 2, gives a list of list... quick fix
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
##' noSpecials(y~us(1|f), delete=FALSE)
##' noSpecials(y~us(1|f), debug=TRUE)
##' noSpecials(y~us+1)  ## should *not* delete unless head of a function
##' noSpecials(~us+1)   ## should work on a one-sided formula!
##' @export
##' @keywords internal
noSpecials <- function(term, delete=TRUE, debug=FALSE) {
    nospec <- noSpecials_(term, delete=delete, debug=debug)
    if (inherits(term, "formula") && length(term) == 3 && is.symbol(nospec)) {
        ## called with two-sided RE-only formula:
        ##    construct response~1 formula
        as.formula(substitute(R~1,list(R=nospec)),
                   env=environment(term))
    } else {
        nospec
    }
}

noSpecials_ <- function(term,delete=TRUE, debug=FALSE) {
    if (debug) print(term)
    if (!anySpecial(term)) return(term)
    if (length(term)==1) return(term)  ## 'naked' specials
    if (isSpecial(term)) {
        if(delete) {
            return(NULL)
        } else { ## careful to return  (1|f) and not  1|f:
            return(substitute((TERM), list(TERM = term[[2]])))
        }
    } else {
        if (debug) print("not special")
        nb2 <- noSpecials_(term[[2]], delete=delete, debug=debug)
        nb3 <- if (length(term)==3) {
                   noSpecials_(term[[3]], delete=delete, debug=debug)
               } else NULL
        if (is.null(nb2)) {
            return(nb3)
        } else if (is.null(nb3)) {
            if (length(term)==2 && identical(term[[1]], quote(`~`))) { ## special case for one-sided formula
                term[[2]] <- nb2
                return(term)
            } else {
                return(nb2)
            }
        } else {  ## neither term completely disappears
            term[[2]] <- nb2
            term[[3]] <- nb3
            return(term)
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
##' extractForm(~offset(x),quote(offset))
##' @export
##' @keywords internal
extractForm <- function(term,value) {
    if (!inForm(term,value)) return(NULL)
    if (is.name(term) || !is.language(term)) return(NULL)
    if (identical(head(term),value)) {
        return(list(term))
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

# n.b. won't work for terms with more than 2 args ...
# @examples
# replaceForm(quote(a(b+x*c(y,z))),quote(y),quote(R))
# ss <- ~(1 | cask:batch) + (1 | batch)
# replaceForm(ss,quote(cask:batch),quote(batch:cask))
replaceForm <- function(term,target,repl) {
    if (identical(term,target)) return(repl)
    if (!inForm(term,target)) return(term)
    if (length(term) == 2) {
        return(substitute(OP(x),list(OP=term[[1]],x=replaceForm(term[[2]],target,repl))))
    }
    return(substitute(OP(x,y),list(OP=term[[1]],
                                   x=replaceForm(term[[2]],target,repl),
                                   y=replaceForm(term[[3]],target,repl))))
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
check_dots <- function(..., action="stop") {
    L <- list(...)
    if (length(L)>0) {
        FUN <- get(action)
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

