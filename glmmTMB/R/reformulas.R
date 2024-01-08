## backward compat (copied from lme4)
if ((getRversion()) < "3.2.1") {
    lengths <- function (x, use.names = TRUE) vapply(x, length, 1L, USE.NAMES = use.names)
}

if (getRversion() < "4.0.0") {
    deparse1 <- function (expr, collapse = " ", width.cutoff = 500L, ...) {
        paste(deparse(expr, width.cutoff, ...), collapse = collapse)
    }
}

#' expand double-bar RE notation by splitting
#' @param term a formula term
#' @rdname formfuns
#' @export
expandDoubleVert <- function(term) {
    frml <- formula(substitute(~x,list(x=term[[2]])))
    ## need term.labels not all.vars to capture interactions too:
    tt <- terms(frml)
    newtrms <- lapply(attr(tt, "term.labels"),
                      function(t) {
                          sumTerms(list(0, toLang(t)))
                      })
    if(attr(tt, "intercept") != 0) {
        newtrms <- c(1, newtrms)
    }
    res <- lapply(newtrms,
           function(t) {
               makeOp(
                   makeOp(t, term[[3]], quote(`|`)),
                   quote(`(`)
               )
           })
    return(res)
}

#' extract right-hand side of a formula
#' @param form a formula object
#' @param as.form (logical) return a formula (TRUE) or as a call/symbolic object (FALSE) ?
#' @export
RHSForm <- function(form, as.form=FALSE) {
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

#' set the right side of a formula
#' @param formula a formula object
#' @param value replacement value for RHS
#' @rdname formfuns
#' @export
`RHSForm<-` <- function(formula,value) {
    formula[[length(formula)]] <- value
    formula
}

#' combine a list of formula terms as a sum
#' @param termList a list of formula terms
#' @rdname formfuns
#' @export
sumTerms <- function(termList) {
    Reduce(function(x,y) makeOp(x,y,op=quote(`+`)),termList)
}

#' extract random effects component of formula
#' @param f a formula
#' @param response include response variable?
#' @param bracket bracket-protect terms?
#' @rdname formfuns
#' @export
reOnly <- function(f, response=FALSE,bracket=TRUE) {
    flen <- length(f)
    f2 <- f[[2]]
    if (bracket)
        f <- lapply(findbars(f), makeOp, quote(`(`)) ## bracket-protect terms
    f <- sumTerms(f)
    if (response && flen==3) {
        form <- makeOp(f2, f, quote(`~`))
    } else {
        form <- makeOp(f, quote(`~`))
    }
    return(form)
}

#' combine unary or binary operator + arguments (sugar for 'substitute')
#' @param x a formula term
#' @param y a formula term (or an operator)
#' @param op an operator
#' @rdname formfuns
## FIXME: would be nice to have multiple dispatch, so
## (arg,op) gave unary, (arg,arg,op) gave binary operator
makeOp <- function(x, y, op=NULL) {
    if (is.null(op) || missing(y)) {  ## unary
        if (is.null(op)) {
            substitute(OP(X),list(X=x,OP=y))
        } else {
            substitute(OP(X),list(X=x,OP=op))
        }
    } else substitute(OP(X,Y), list(X=x,OP=op,Y=y))
}

#' combines the right-hand sides of two formulas, or a formula and a symbol
#' @param f1 formula #1
#' @param f2 formula #2
#' @rdname formfuns
#' @export
#' @examples
#' addForm0(y~x,~1)
#' addForm0(~x,~y)
addForm0 <- function(f1,f2) {
    tilde <- as.symbol("~")
    if (!identical(head(f2),tilde)) {
        f2 <- makeOp(f2,tilde)
    }
    if (length(f2)==3) warning("discarding LHS of second argument")
    RHSForm(f1) <- makeOp(RHSForm(f1),RHSForm(f2),quote(`+`))
    return(f1)
}

#' Combine right-hand sides of an arbitrary number of formulas
#' @param ... arguments to pass through to \code{addForm0}
#' @rdname formfuns
#' @export
addForm <- function(...) {
  Reduce(addForm0,list(...))
}

#' list of specials -- taken from enum.R
findReTrmClasses <- function() {
    c(names(.valid_covstruct), "s")
}

toLang <- function(x) parse(text=x)[[1]]

#' apply
#' @param f a language object (an atom of a formula)
#' expandGrpVar(quote(x*y))
#' expandGrpVar(quote(x/y))
expandGrpVar <- function(f) {
    form <- as.formula(makeOp(f,quote(`~`)))
    mm <- terms(form)
    tl <- attr(mm,"term.labels")
    ## reverse order: f/g -> f + g:f (for lme4/back-compatibility)
    switch_order <- function(x) paste(rev(unlist(strsplit(x, ":"))), collapse = ":")
    if (inForm(f, quote(`/`))) {
        ## vapply adds names; remove them, and reverse order of sub-terms, for back-compatibility ...
        tl <- unname(vapply(tl, switch_order, character(1)))
        tl <- rev(tl)
    }
    res <- lapply(tl, toLang)
    return(res)
}

##' expand interactions/combinations of grouping variables
##'
##' Modeled after lme4:::expandSlash, by Doug Bates. However,
##' all formula operators that apply to factors (\code{*}, \code{/}, \code{+})
##' are applicable: the results are expanded into a list of independent (additive)
##' random effect terms
##' @param bb a list of naked grouping variables, i.e. 1 | f
##' @examples
##' ff <- lme4::findbars(y~1+(x|f/g))
##' expandAllGrpVar(ff)
##' expandAllGrpVar(quote(1|(f/g)/h))
##' expandAllGrpVar(quote(1|f/g/h))
##' expandAllGrpVar(quote(1|f*g))
##' expandAllGrpVar(quote(1|f+g))
##' expandAllGrpVar(quote(a+b|f+g+h*i))
##' ## wish list ... this should be (1|a) + (1|a:b) + (1|a:b:c) + (1|a:b:d) ...
##' ## expandAllGrpVar(quote(a/b/(c+d)))
##' expandAllGrpVar(quote(s(log(d), k = 4)))
##' expandAllGrpVar(quote(s(log(d+1))))
##' @importFrom utils head
##' @rdname formfuns
##' @export
expandAllGrpVar <- function(bb) {
    ## Return the list of expanded terms (/, *, ?)
    if (!is.list(bb))
        expandAllGrpVar(list(bb))
    else {
        for (i in seq_along(bb)) {
            esfun <- function(x) {
                if (length(x)==1 || !anySpecial(x, "|")) return(x)
                if (length(x)==2) {
                        ## unary operator such as diag(1|f/g)
                        ## return diag(...) + diag(...) + ...
                        return(lapply(esfun(x[[2]]),
                                      makeOp, y=head(x)))
                }
                if (length(x)==3) {
                    ## binary operator
                    if (x[[1]]==quote(`|`)) {
                        return(lapply(expandGrpVar(x[[3]]),
                                      makeOp, x=x[[2]], op=quote(`|`)))
                    } else {
                        return(x)
                        ## return(x) would be nice, but in that case x gets evaluated
                        ## return(setNames(makeOp(esfun(x[[2]]), esfun(x[[3]]),
                        ##  op=x[[1]]), names(x)))
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
head.name <- function(x, ...) { x }

## TEST: does this work as a drop-in replacement for lme4::findbars
## if default.special = NULL?
## (would replace current expandDoubleVerts machinery)

##' Find and process random effects terms
##'
##' @param term a formula or piece of a formula
##' @param debug (logical) debug?
##' @param specials list of special terms
##' @param default.special character: special to use for parenthesized terms - i.e. random effects terms with unspecified structure
##' @param expand_doublevert_method method for handling \code{||} operator: split into separate terms or replace by \code{diag}? Inherited from \emph{previous call where it was specified}.
##' 1. atom (not a call or an expression): NULL
##' 2. special, i.e. foo(...) where "foo" is in specials: return term
##' 3. parenthesized term: \emph{if} the head of the head is | (i.e.
##'    it is of the form (xx|gg), then convert it to the default
##'    special type; we won't allow pathological cases like
##'    ((xx|gg)) ... [can we detect them?]
##' @examples
##' splitForm(quote(us(x,n=2)))
##' findbars_x(~ 1 + (x + y || g), expand_doublevert_method = "diag_special")
##' findbars_x(~ 1 + (x + y || g), expand_doublevert_method = "split")
##' findbars_x(~ 1 + (1 | f) + (1 | g))
##' findbars_x(~ 1 + (1 | f) + (1 | g))
##' findbars_x(~ 1 + (1|h) + (x + y || g), expand_doublevert_method = "split")
##' findbars_x(~ 1 + (1|Subject))
##' findbars_x(~ (1||Subject))
##' findbars_x(~ (1|Subject))
##' findbars_x(~ (1|Subject), default.special = NULL)
##' findbars_x(~ 1 + x)
##' findbars_x(~ s(x, bs = "tp"))
##' findbars_x(y ~ a + log(b) + s(x, bs = "tp") + s(y, bs = "gp"),
##'    target = "s", default.special = NULL)
##' @rdname formfuns
##' @export
findbars_x <- function(term,
                debug=FALSE,
                specials=character(0),
                default.special="us",
                target = '|',
                expand_doublevert_method = c("diag_special", "split")) {

    expand_doublevert_method <- match.arg(expand_doublevert_method)
    
    ds <- if (is.null(default.special)) {
              NULL
          } else {
              ## convert default special char to symbol (less ugly way?)
              eval(substitute(as.name(foo),list(foo=default.special)))
          }

    ## base function
    ## defining internally in this way makes debugging slightly
    ## harder, but (1) allows easy propagation of the top-level
    ## arguments down the recursive chain; (2) allows the top-level
    ## expandAllGrpVar() operation (which also handles cases where
    ## a naked term rather than a list is returned)
    
    fbx <- function(term) {
        if (is.name(term) || !is.language(term)) return(NULL)
        if (list(term[[1]]) %in% lapply(specials,as.name)) {
            if (debug) cat("special: ",deparse(term),"\n")
            return(term)
        }
        if (head(term) == as.name(target)) {  ## found x | g
            if (debug) {
                tt <- if (target == '|') "bar" else sprintf('"%s"', target)
                cat(sprintf("%s term: %s\n", tt, deparse(term)))
            }
            if (is.null(ds)) return(term)
            return(makeOp(term, ds))
        }
        if (head(term) == as.name("||")) {
            if (expand_doublevert_method == "diag_special") {
                return(makeOp(makeOp(term[[2]], term[[3]],
                                     op = quote(`|`)),
                              as.name("diag")))
            }
            if (expand_doublevert_method == "split") {
                ## need to return *multiple* elements
                return(lapply(expandDoubleVert(term), fbx))
            } 
            stop("unknown doublevert method ", expand_doublevert_method)
        }
        if (head(term) == as.name("(")) {  ## found (...)
            if (debug) cat("paren term:",deparse(term),"\n")
            return(fbx(term[[2]]))
        }
        stopifnot(is.call(term))
        if (length(term) == 2) {
            ## unary operator, decompose argument
            if (debug) cat("unary operator:",deparse(term[[2]]),"\n")
            return(fbx(term[[2]]))
        }
        ## binary operator, decompose both arguments
        f2 <- fbx(term[[2]])
        f3 <- fbx(term[[3]])

        if (debug) { cat("binary operator:",deparse(term[[2]]),",",
                         deparse(term[[3]]),"\n")
                         cat("term 2: ", deparse(f2), "\n")
                         cat("term 3: ", deparse(f3), "\n")
        }
        c(f2, f3)
    }

    fbx_term <- fbx(term)
    if (debug) cat("fbx(term): ", deparse(fbx_term))
    expandAllGrpVar(fbx_term)

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
##' splitForm(~1+s(x, bs = "tp"))
##'
##' @author Steve Walker
##' @importFrom lme4 nobars
##' @export
splitForm <- function(formula,
                      defaultTerm="us",
                      allowFixedOnly=TRUE,
                      allowNoSpecials=TRUE,
                      debug=FALSE,
                      specials = findReTrmClasses()) {

    ## logic:

    ## string for error message *if* specials not allowed
    ## (probably package-specific)
    noSpecialsAlt <- "lmer or glmer"

    ## formula <- expandDoubleVerts(formula)
    ## split formula into separate
    ## random effects terms
    ## (including special terms)

    fbxx <- findbars_x(formula, debug, specials)
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

    ## nobars() will get rid of any *naked* RE terms
    ## FIXME ... let noSpecials handle naked bar-terms if desired ?
    ## (would adding "|" to reTrmClasses work?)
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
##' noSpecials(~us(1|f)+1)   ## should work on a one-sided formula!
##' noSpecials(~s(stuff) + a + b, specials = "s")
##' @export
##' @keywords internal
noSpecials <- function(term, delete=TRUE, debug=FALSE, specials = findReTrmClasses()) {
    nospec <- noSpecials_(term, delete=delete, debug=debug, specials = specials)
    if (inherits(term, "formula") && length(term) == 3 && is.symbol(nospec)) {
        ## called with two-sided RE-only formula:
        ##    construct response~1 formula
        as.formula(substitute(R~1,list(R=nospec)),
                   env=environment(term))
        ## FIXME::better 'nothing left' handling
    } else if (is.null(nospec)) {
        ~1
    } else {
        nospec
    }
}

noSpecials_ <- function(term, delete=TRUE, debug=FALSE, specials = findReTrmClasses()) {
    if (debug) print(term)
    if (!anySpecial(term, specials)) return(term)
    if (length(term)==1) return(term)  ## 'naked' specials
    if (isSpecial(term, specials)) {
        if(delete) {
            return(NULL)
        } else { ## careful to return  (1|f) and not  1|f:
            return(substitute((TERM), list(TERM = term[[2]])))
        }
    } else {
        if (debug) print("not special")
        nb2 <- noSpecials_(term[[2]], delete=delete, debug=debug, specials = specials)
        nb3 <- if (length(term)==3) {
                   noSpecials_(term[[3]], delete=delete, debug=debug, specials = specials)
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

isSpecial <- function(term, specials = findReTrmClasses()) {
    if(is.call(term)) {
        ## %in% doesn't work (requires vector args)
        for(cls in specials) {
            if(term[[1]] == cls) return(TRUE)
        }
    }
    FALSE
}

isAnyArgSpecial <- function(term, specials = findReTrmClasses()) {
    for(tt in term)
        if(isSpecial(tt, specials)) return(TRUE)
    FALSE
}

## This could be in principle be fooled by a term with a matching name
## but this case is caught in noSpecials_() where we test for length>1
anySpecial <- function(term, specials=findReTrmClasses()) {
    any(specials %in% all.names(term))
}

##' test whether a formula contains a particular element?
##' @rdname formfuns
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
inForm <- function(form, value) {
    if (any(sapply(form,identical,value))) return(TRUE)
    if (all(sapply(form,length)==1)) return(FALSE)
    return(any(vapply(form,inForm,value,FUN.VALUE=logical(1))))
}

##' extract terms with a given head from an expression/formula
##' @rdname formfuns
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
##' @rdname formfuns
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


##' drop terms matching a particular value from an expression
##' @rdname formfuns
## from Gabor Grothendieck: recursive solution
## http://stackoverflow.com/questions/40308944/removing-offset-terms-from-a-formula
##' @param x formula
##' @param value term to remove from formula
##' @param preserve (integer) retain the specified occurrence of "value"
##' @examples
##' drop.special(x~a + b+ offset(z))
##' @export
##' @keywords internal
drop.special <- function(x, value=quote(offset), preserve = NULL) {
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

#' replace a component of a call/parse tree
#' n.b. won't work for terms with more than 2 args ...
#' @rdname formfuns
#' @export
#' @examples
#' replaceForm(quote(a(b+x*c(y,z))),quote(y),quote(R))
#' ss <- ~(1 | cask:batch) + (1 | batch)
#' replaceForm(ss,quote(cask:batch),quote(batch:cask))
#' replaceForm(ss, quote(`:`), quote(`%:%`))
replaceForm <- function(term,target,repl) {
    if (identical(term,target)) return(repl)
    if (!inForm(term,target)) return(term)
    if (length(term) == 2) {
        return(substitute(OP(x),list(OP=replaceForm(term[[1]],target,repl),
                                     x=replaceForm(term[[2]],target,repl))))
    }
    return(substitute(OP(x,y),list(OP=replaceForm(term[[1]],target,repl),
                                   x=replaceForm(term[[2]],target,repl),
                                   y=replaceForm(term[[3]],target,repl))))
}

##' @noRd
##' @examples
##' no_specials(y ~ 1 + s(x) + (f|g))
no_specials <- function(term, specials = c("|", "||", "s")) {
    if (is.list(term)) {
        return(lapply(term, no_specials))
    }
    for (ss in specials) {
        if (identical(head(term), as.name(ss))) return(term)
    }
    if (length(term) == 3) stop("don't know what to do")
    return(no_specials(term[[2]], specials))
}


##' Substitute safe chars (+) for specials (for use in \code{model.frame})
##' (Generalized from \code{lme4}'s \code{subbars} function.)
##' @param term formula or term in a formula
##' @param specials names of specials to process
##' @param keep_args number of arguments to retain (matching \code{specials})
##' @return a term or formula with specials replaced by \code{+} (and extra arguments dropped)
##' @keywords internal
##' @examples
##' sub_specials( ~ s(a, k=4))
##' sub_specials( ~ (1|x) + (a + b || y) + s(a, k=4))
##' sub_specials(Reaction ~ s(Days) + (1 + Subject))
##' sub_specials(~ s(cos((y^2*3)/2), bs = "tp"))
##' @export
sub_specials <- function (term,
                          specials = c("|", "||", "s"),
                          keep_args = c(2L, 2L, NA_integer_)) {
    if (is.name(term) || !is.language(term)) 
        return(term)
    ## previous version recursed immediately for unary operators,
    ## (we were only interested in `|`(x,y) and `||`(x,y))
    ## but here s(x) needs to be processed ...
    for (i in seq_along(specials)) {
        if (is.call(term) && term[[1]] == as.name(specials[i])) {
            if (is.na(keep_args[i])) {
                ## keep only *unnamed* args
                if (!is.null(names(term))) {
                    term <- term[names(term)==""]
                }
            } else {
                term <- term[1:(1+keep_args[i])]
            }
            term[[1]] <- as.name("+")
            ## converts s(x) to +x, which is ugly, but
            ##  formula can handle repeated '+'
            ## discard additional arguments (e.g for s(x, ...))
            ## (fragile re: order??)
        }
    }
    for (j in 2:length(term)) {
        term[[j]] <- sub_specials(term[[j]],
                                  specials = specials,
                                  keep_args = keep_args)
    }
    term
}
