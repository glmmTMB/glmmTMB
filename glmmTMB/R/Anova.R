## Type II and III tests for linear, generalized linear, and other models (J. Fox)
## most of what's below is copied from car::Anova.R
## main changes are (1) absence of F-test (K-R, Satterthwaite df) capability;
## (2) use of [[component]] to pick out relevant fixed effect parameters/v-cov matrix

## copied unchanged (?); unexported utilities from car
responseName.default <- function (model, ...) deparse(attr(terms(model), "variables")[[2]])

term.names.default <- function (model, component="cond", ...) {
    term.names <- labels(terms(model, component=component))
    if (has.intercept(model)) c("(Intercept)", term.names)
    else term.names
}

has.intercept <- function (model, ...) {
	UseMethod("has.intercept")
}

ConjComp <- function(X, Z = diag( nrow(X)), ip = diag(nrow(X))) {
    ## This function by Georges Monette
    ## finds the conjugate complement of the proj of X in span(Z) wrt
    ##    inner product ip
    ## - assumes Z is of full column rank
    ## - projects X conjugately wrt ip into span Z
    xq <- qr(t(Z) %*% ip %*% X)
    if (xq$rank == 0) return(Z)
    Z %*% qr.Q(xq, complete = TRUE) [ ,-(1:xq$rank)] 
}

relatives <- function(term, names, factors){
  is.relative <- function(term1, term2) {
    all(!(factors[,term1]&(!factors[,term2])))
  }
  if(length(names) == 1) return(NULL)
  which.term <- which(term==names)
  (1:length(names))[-which.term][sapply(names[-which.term], 
                                        function(term2) is.relative(term, term2))]
}

## modified
has.intercept.glmmTMB <- function (model, component="cond", ...) {
    nms <- names(fixef(model)[[component]])
    any(grepl("\\(Intercept\\)",nms))
}

##' @rdname downstream_methods
##' @rawNamespace if(getRversion() >= "3.6.0") {
##'      S3method(car::Anova, glmmTMB)
##' } else {
##'   export(Anova.glmmTMB)
##' }

##' @param vcov. variance-covariance matrix (usually extracted automatically)
##' @param test.statistic unused: only valid choice is "Chisq" (i.e., Wald chi-squared test)
##' @param singular.ok OK to do ANOVA with singular models (unused) ?
##' @param type  type of test, \code{"II"}, \code{"III"}, \code{2}, or \code{3}.  Roman numerals are equivalent to the corresponding Arabic numerals. See \code{\link[car]{Anova}} for details.

Anova.glmmTMB <- function (mod, type = c("II", "III", 2, 3),
                           test.statistic = c("Chisq","F"),
                           component="cond",
                           vcov. = vcov(mod)[[component]], singular.ok, ...) 
{

    ff <- fixef(mod)[[component]]
    if (trivialFixef(names(ff),component)) {
        stop(sprintf("trivial fixed effect for component %s: can't compute Anova table", sQuote(component)))
    }
    test.statistic <- match.arg(test.statistic)
    if (test.statistic=="F") {
        stop("F tests currently unavailable")
    }
    if (is.function(vcov.)) 
        vcov. <- vcov.(mod)
    type <- as.character(type)
    type <- match.arg(type)
    if (missing(singular.ok)) 
        singular.ok <- type == "2" || type == "II"
    afun <- switch(type,
                   `2` = , II = Anova.II.glmmTMB,
                   `3` = , III = Anova.III.glmmTMB)
    afun(mod, vcov., test=test.statistic, singular.ok = singular.ok,
         component = component)
}

## defined as a function, not a method, so we can hand the object
## off to car::linearHypothesis.default (not exported)

linearHypothesis_glmmTMB <- function (model, hypothesis.matrix,
                      rhs = NULL, test = c("Chisq", "F"),
                      vcov. = NULL, singular.ok = FALSE, verbose = FALSE, 
                      coef. = NULL, component="cond", ...)
{
    ## what's the least ugly way to do this?
    ## match.call?
    test <- match.arg(test)
    ## call linearHypothesis.default (not exported)
    if (!requireNamespace("car")) {
        stop("please install (if necessary) and load the car package")
    }
    if (utils::packageVersion("car")<"3.0.6") {
        stop("please install a more recent version of the car package (>= 3.0.6)")
    }
    car::linearHypothesis(model=model,
             hypothesis.matrix=hypothesis.matrix,
             rhs=rhs,
             test=test,
             vcov. = vcov(model)[[component]],
             singular.ok = FALSE,
             verbose = verbose,
             coef. = fixef(model)[[component]],
             ...)
}                  
    
Anova.II.glmmTMB <- function(mod, vcov., singular.ok=TRUE, test="Chisq",
                             component="cond", ...){

    ## would feel cleaner to have this external, but it uses
    ##  lots of variable from the function environment ...
    hyp.term <- function(term) {
        which.term <- which(term==names)
        subs.term <- which(assign==which.term)
        relatives <- relatives(term, names, fac)
        subs.relatives <- NULL
        for (relative in relatives) 
            subs.relatives <- c(subs.relatives, which(assign==relative))
        hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
        hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
        hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
        hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]       
        hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) {
                               hyp.matrix.2
                           } else {
                               t(ConjComp(t(hyp.matrix.1),
                                          t(hyp.matrix.2), vcov.))
                           }
        hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
                                                  function(x) all(x == 0)), , drop=FALSE]
        if (nrow(hyp.matrix.term) == 0)
            return(c(statistic=NA, df=0))            
        hyp <- linearHypothesis_glmmTMB(mod, hyp.matrix.term, 
                                        vcov.=vcov.,
                                        singular.ok=singular.ok,
                                        test=test,
                                        component=component, ...)
        if (test == "Chisq") return(c(statistic=hyp$Chisq[2], df=hyp$Df[2]))
        else return(c(statistic=hyp$F[2], df=hyp$Df[2], res.df=hyp$Res.Df[2]))
    } ## hyp.term()

    ## may be irrelevant, glmmTMB doesn't currently handle aliased terms?
    not.aliased <- !is.na(fixef(mod)[[component]])
    if (!singular.ok && !all(not.aliased))
        stop("there are aliased coefficients in the model")
    fac <- attr(terms(mod, component=component), "factors")
    intercept <- has.intercept(mod)
    p <- length(fixef(mod)[[component]])
    I.p <- diag(p)
    if (!missing(vcov.)){
        vcov. <- vcov(mod, complete=FALSE)[[component]]
    }
    assign <- attr(model.matrix(mod, component=component), "assign")
    assign[!not.aliased] <- NA
    names <- term.names.default(mod, component=component)
    if (intercept) names <- names[-1]
    n.terms <- length(names)
    p <- teststat <- df <- res.df <- rep(0, n.terms)
    for (i in seq_len(n.terms)) {
        hyp <- hyp.term(names[i])
        teststat[i] <- abs(hyp["statistic"])
        df[i] <- abs(hyp["df"])
        res.df[i] <- hyp["res.df"]
        p[i] <- pchisq(teststat[i], df[i], lower.tail=FALSE) 
    } 
    result <- data.frame(teststat, df, p)
    row.names(result) <- names
    names(result) <- c ("Chisq", "Df", "Pr(>Chisq)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c("Analysis of Deviance Table (Type II Wald chisquare tests)\n", 
                                 paste("Response:", responseName.default(mod)))
    return(result)
}

Anova.III.glmmTMB <- function(mod, vcov., singular.ok=FALSE, test="Chisq",
                              component="cond", ...){
    intercept <- has.intercept(mod)
    p <- length(fixef(mod)[[component]])
    I.p <- diag(p)
    names <- term.names.default(mod, component=component)
    n.terms <- length(names)
    assign <- attr(model.matrix(mod, component=component), "assign")
    p <- teststat <- df <- res.df <- rep(0, n.terms)
    if (intercept) df[1] <- 1
    not.aliased <- !is.na(fixef(mod)[[component]])
    if (!singular.ok && !all(not.aliased))
        stop("there are aliased coefficients in the model")
    if (!missing(vcov.)){
        vcov. <- vcov(mod, complete=FALSE)[[component]]
    }
    for (term in seq_len(n.terms)){
        subs <- which(assign == term - intercept)
        hyp.matrix <- I.p[subs,,drop=FALSE]
        hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
        hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]        
        if (nrow(hyp.matrix) == 0){
            teststat[term] <- NA
            df[term] <- 0
            p[term] <- NA
        }
        else {
            hyp <- linearHypothesis_glmmTMB(mod, hyp.matrix, test=test,
                                            vcov.=vcov., singular.ok=singular.ok,
                                            component=component, ...)
            teststat[term] <-  hyp$Chisq[2] 
            df[term] <- abs(hyp$Df[2])
            p[term] <- pchisq(teststat[term], df[term], lower.tail=FALSE) 
        }
        result <- data.frame(teststat, df, p)
        row.names(result) <- names
        names(result) <- c ("Chisq", "Df", "Pr(>Chisq)")
        class(result) <- c("anova", "data.frame")
        attr(result, "heading") <- c("Analysis of Deviance Table (Type III Wald chisquare tests)\n", 
                                     paste("Response:", responseName.default(mod)))
    }
    result
}
