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

## roxygen complains if this is not exported ... ?
## modified
#' @export
has.intercept.glmmTMB <- function (model, component="cond", ...) {
    nms <- names(fixef(model)[[component]])
    any(grepl("\\(Intercept\\)",nms))
}

## n.b. rawNamespace spec must now be on a single line!

##' @rdname downstream_methods
##' @rawNamespace if(getRversion() >= "3.6.0") { S3method(car::Anova, glmmTMB) } else { export(Anova.glmmTMB) }
##' @param vcov. variance-covariance matrix (usually extracted automatically); not
##' currently supported together with \code{ddf != "asymptotic"}
##' @param test.statistic \code{"Chisq"} (default; a Wald chi-squared test) or \code{"F"}
##' (only available together with \code{ddf != "asymptotic"}, see \code{ddf} below). An explicit
##' \code{test.statistic = "Chisq"} combined with \code{ddf != "asymptotic"} is an error, since
##' Kenward-Roger/Satterthwaite always produce an F table
##' @param singular.ok OK to do ANOVA with singular models (unused) ?
##' @param type  type of test, \code{"II"}, \code{"III"}, \code{2}, or \code{3}.  Roman numerals are equivalent to the corresponding Arabic numerals. See \code{\link[car]{Anova}} for details.
##' @param include.rankdef.cols include all columns of a rank-deficient model matrix?
##' @param ddf denominator degrees-of-freedom calculation, as in \code{\link{summary.glmmTMB}}
##' and \code{\link{anova.glmmTMB}}. The default \code{"asymptotic"} gives the classical Wald
##' chi-squared table; \code{"kenward-roger"} or \code{"satterthwaite"} instead give an F-ratio
##' table, with each term's denominator df computed via the Kenward-Roger or Satterthwaite
##' approximation (see \code{\link{dof_KR}}, \code{\link{dof_satt}}). \code{"kenward-roger"}
##' additionally requires a family with an estimated dispersion parameter, and throws an error
##' for families such as \code{binomial} or \code{poisson} that lack one (\code{"satterthwaite"}
##' has no such restriction). Not currently supported together with a user-supplied \code{vcov.},
##' \code{component != "cond"}, or models with aliased/rank-deficient or \code{map}-fixed
##' conditional coefficients.

Anova.glmmTMB <- function (mod, type = c("II", "III", 2, 3),
                           test.statistic = c("Chisq","F"),
                           component="cond",
                           vcov. = vcov(mod)[[component]],
                           singular.ok,
                           include.rankdef.cols = FALSE,
                           ddf = c("asymptotic", "kenward-roger", "satterthwaite"),
                           ...) {

    ff <- fixef(mod)[[component]]
    if (trivialFixef(names(ff),component)) {
        stop(sprintf("trivial fixed effect for component %s: can't compute Anova table", sQuote(component)))
    }
    ddf <- match.arg(ddf)
    user_test_statistic <- !missing(test.statistic)
    test.statistic <- match.arg(test.statistic)
    if (test.statistic == "F" && ddf == "asymptotic") {
        stop("F tests require ddf='kenward-roger' or ddf='satterthwaite' ",
             "(F tests are not available for ddf='asymptotic')")
    }
    user_vcov <- !missing(vcov.)
    if (ddf != "asymptotic") {
        if (component != "cond") {
            stop("ddf is currently only supported for component = 'cond'")
        }
        if (user_vcov) {
            stop("a user-supplied 'vcov.' cannot be combined with ddf != 'asymptotic': ",
                 "Kenward-Roger/Satterthwaite need the model's own REML/ML ",
                 "variance-parameter uncertainty, not an arbitrary covariance matrix")
        }
        ## ddf != "asymptotic" always produces an F table (see
        ## Anova.II/III.glmmTMB's ddf branch, which ignores 'test'
        ## entirely); an explicit test.statistic="Chisq" would silently
        ## be overridden, so reject that combination instead, and quietly
        ## default to "F" otherwise
        if (user_test_statistic && test.statistic == "Chisq") {
            stop("test.statistic='Chisq' cannot be combined with ddf != 'asymptotic': ",
                 "Kenward-Roger/Satterthwaite always produce F tests; ",
                 "omit test.statistic or set it to 'F'")
        }
        test.statistic <- "F"
        check_ddf(mod, ddf)
    }
    if (is.function(vcov.))
        vcov. <- vcov.(mod)
    ## coefficients fixed via 'map' are known constants (zero variance);
    ## adjust the default vcov so it aligns with the full coefficient
    ## vector -- a user-supplied matrix is trusted as-is
    if (!user_vcov) vcov. <- pad_mapped_vcov(mod, vcov., component)
    type <- as.character(type)
    type <- match.arg(type)
    if (missing(singular.ok))
        singular.ok <- type == "2" || type == "II"
    afun <- switch(type,
                   `2` = , II = Anova.II.glmmTMB,
                   `3` = , III = Anova.III.glmmTMB)
    afun(mod, vcov., test=test.statistic, singular.ok = singular.ok,
         component = component, include.rankdef.cols = include.rankdef.cols,
         ddf = ddf)
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
             vcov. = vcov.,
             singular.ok = singular.ok,
             verbose = verbose,
             coef. = fixef(model)[[component]],
             ...)
}                  
    
## shared setup for the ddf-based (Kenward-Roger/Satterthwaite) F-test path
## in Anova.II/III.glmmTMB: adds the Anova-specific "no aliased/map-fixed
## coefficients" guard (check_ddf(), called once already by Anova.glmmTMB(),
## validates everything else -- REML, family, presence of random effects)
## and precomputes the potentially expensive Kenward-Roger adjusted vcov
## once so it isn't recomputed for every term
.Anova_ddf_setup <- function(mod, component, ddf, not.aliased) {
    if (ddf == "asymptotic") return(NULL)
    if (!all(not.aliased)) {
        stop("ddf='", ddf, "' is not currently supported for models with ",
             "aliased/rank-deficient conditional coefficients; use ddf='asymptotic'")
    }
    unadjusted_vcov <- stats::vcov(mod)[[component]]
    ## map-fixed coefficients keep their fixed *value* in fixef() (unlike
    ## aliased/rank-deficient ones, which show up as NA there and are
    ## already caught above by not.aliased) -- in the raw, unpadded vcov()
    ## used here they instead show up as NA variance (pad_mapped_vcov(),
    ## used by the asymptotic path, turns that NA into an explicit 0; this
    ## path bypasses pad_mapped_vcov() since a padded/zero-variance matrix
    ## isn't meaningful input for Kenward-Roger/Satterthwaite). Left
    ## unguarded, this reaches .vcov_kenward_adjusted()'s linear algebra
    ## and fails with an opaque eigen() error instead
    if (anyNA(diag(unadjusted_vcov))) {
        stop("ddf='", ddf, "' is not currently supported for models with ",
             "map-fixed conditional coefficients; use ddf='asymptotic'")
    }
    has_random <- hasRandom(mod)
    list(hasRandom = has_random,
         beta = fixef(mod)[[component]],
         unadjusted_vcov = unadjusted_vcov,
         residual_df = stats::df.residual(mod),
         adjusted_vcov = if (ddf == "kenward-roger" && has_random) .vcov_kenward_adjusted(mod) else NULL)
}

## per-term F-test given an already-built hypothesis matrix; used in place
## of linearHypothesis_glmmTMB() when ddf != "asymptotic". Dispatches to the
## same joint Kenward-Roger/Satterthwaite/no-random-effects machinery
## anova.glmmTMB() uses for model-comparison F-tests (denom_df.R), just
## given a hypothesis matrix directly instead of a pair of nested models
.Anova_ddf_test <- function(mod, hyp.matrix, ddf, ddf_info) {
    if (nrow(hyp.matrix) == 0) {
        return(c(Fstat = NA_real_, ndf = 0, ddf = NA_real_, pval = NA_real_))
    }
    res <- if (!ddf_info$hasRandom) {
        .wald_joint_test(ddf_info$unadjusted_vcov, hyp.matrix, ddf_info$beta,
                         ddf = ddf_info$residual_df)
    } else if (ddf == "kenward-roger") {
        .KR_adjust_joint(ddf_info$adjusted_vcov, ddf_info$unadjusted_vcov,
                         hyp.matrix, ddf_info$beta)
    } else {
        .satt_adjust_joint(mod, hyp.matrix)
    }
    c(Fstat = res$Fstat, ndf = res$ndf, ddf = res$ddf, pval = res$p.value)
}

## rbind() a list of per-term result rows, with the right column names even
## when there are zero terms (e.g. an intercept-only zi/disp formula) --
## do.call(rbind, list()) returns NULL rather than a 0-row matrix, which
## would otherwise break the rows[, "statistic"]-style indexing downstream
.rbind_terms <- function(row_list, cols) {
    if (length(row_list) == 0) {
        return(matrix(numeric(0), nrow = 0, ncol = length(cols),
                      dimnames = list(NULL, cols)))
    }
    do.call(rbind, row_list)
}

## assemble the Type II/III F-ratio result table for ddf != "asymptotic",
## given a matrix of Fstat/ndf/ddf/pval rows (one per term); shared between
## Anova.II.glmmTMB and Anova.III.glmmTMB since this code path is new (no
## legacy row/column conventions from the original Wald chi-squared path to
## preserve, unlike the two asymptotic-path table-builders below)
.Anova_ddf_result_table <- function(names, ddf, rows, type, response) {
    result <- as.data.frame(rows[, c("Fstat", "ndf", "ddf", "pval"), drop = FALSE])
    names(result) <- c("F", "Num Df", "Den Df", "Pr(>F)")
    row.names(result) <- names
    method <- if (ddf == "kenward-roger") "Kenward-Roger" else "Satterthwaite"
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c(
        sprintf("Analysis of Deviance Table (Type %s %s F tests)\n", type, method),
        paste("Response:", response))
    result
}

Anova.II.glmmTMB <- function(mod, vcov., singular.ok=TRUE, test="Chisq",
                             component="cond", include.rankdef.cols = FALSE,
                             ddf = "asymptotic", ...){

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
        ## hypothesis rows involving only map-fixed coefficients (known
        ## constants, zero variance) are untestable; drop them so the
        ## term gets an NA row instead of a singular-matrix error
        ## (guard against NA variances from user-supplied vcov, e.g.
        ## include_nonest=TRUE, which would make any(zv) return NA)
        zv <- !is.na(diag(vcov.)) & diag(vcov.) == 0
        if (any(zv) && nrow(hyp.matrix.term) > 0) {
            hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1,
                                          function(x) all(x[!zv] == 0)), ,
                                          drop=FALSE]
        }
        if (ddf != "asymptotic") {
            return(.Anova_ddf_test(mod, hyp.matrix.term, ddf, ddf_info))
        }
        ## 'test' can only be "Chisq" here: Anova.glmmTMB() rejects
        ## test.statistic="F" unless ddf != "asymptotic", which is handled above
        if (nrow(hyp.matrix.term) == 0)
            return(c(statistic=NA, df=0))
        hyp <- linearHypothesis_glmmTMB(mod, hyp.matrix.term,
                                        vcov.=vcov.,
                                        singular.ok=singular.ok,
                                        test=test,
                                        component=component, ...)
        return(c(statistic=hyp$Chisq[2], df=hyp$Df[2]))
    } ## hyp.term()

    not.aliased <- !is.na(fixef(mod)[[component]])
    if (!singular.ok && !all(not.aliased))
        stop("there are aliased coefficients in the model")
    fac <- attr(terms(mod, component=component), "factors")
    intercept <- has.intercept(mod)
    p <- length(fixef(mod)[[component]])
    I.p <- diag(p)
    ## FIXME:: missing or !missing ???
    if (missing(vcov.)){
        vcov. <- vcov(mod, complete=FALSE)[[component]]
    }
    vcov. <- vcov.[not.aliased, not.aliased]
    assign <- attr(model.matrix(mod, component=component, include_rankdef = include.rankdef.cols), "assign")
    assign[!not.aliased] <- NA
    names <- term.names.default(mod, component=component)
    if (intercept) names <- names[-1]
    n.terms <- length(names)
    ddf_info <- .Anova_ddf_setup(mod, component, ddf, not.aliased)
    rows <- .rbind_terms(lapply(names, hyp.term),
                         if (ddf == "asymptotic") c("statistic", "df")
                         else c("Fstat", "ndf", "ddf", "pval"))
    if (ddf != "asymptotic") {
        return(.Anova_ddf_result_table(names, ddf, rows, type = "II",
                                       response = responseName.default(mod)))
    }
    teststat <- abs(rows[, "statistic"])
    df <- abs(rows[, "df"])
    result <- data.frame(teststat, df, p = pchisq(teststat, df, lower.tail = FALSE))
    row.names(result) <- names
    names(result) <- c("Chisq", "Df", "Pr(>Chisq)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c("Analysis of Deviance Table (Type II Wald chisquare tests)\n",
                                 paste("Response:", responseName.default(mod)))
    return(result)
}

Anova.III.glmmTMB <- function(mod, vcov., singular.ok=FALSE, test="Chisq",
                              component="cond", include.rankdef.cols = FALSE,
                              ddf = "asymptotic", ...){
    intercept <- has.intercept(mod)
    p <- length(fixef(mod)[[component]])
    I.p <- diag(p)
    names <- term.names.default(mod, component=component)
    n.terms <- length(names)
    not.aliased <- !is.na(fixef(mod)[[component]])
    if (!singular.ok && !all(not.aliased))
        stop("there are aliased coefficients in the model")
    if (missing(vcov.)){
        vcov. <- vcov(mod, complete=FALSE)[[component]]
    }
    vcov. <- vcov.[not.aliased, not.aliased]
    assign <- attr(model.matrix(mod, component=component, include_rankdef = include.rankdef.cols), "assign")
    assign[!not.aliased] <- NA

    term.hyp.matrix <- function(term) {
        subs <- which(assign == term - intercept)
        hyp.matrix <- I.p[subs,,drop=FALSE]
        hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
        hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]
        ## hypothesis rows involving only map-fixed coefficients (known
        ## constants, zero variance -- e.g. the ordinal family's intercept)
        ## are untestable; drop them so the term gets an NA row instead of a
        ## singular-matrix error (guard against NA variances from
        ## user-supplied vcov, e.g. include_nonest=TRUE, which would make
        ## any(zv) return NA)
        zv <- !is.na(diag(vcov.)) & diag(vcov.) == 0
        if (any(zv) && nrow(hyp.matrix) > 0) {
            hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1,
                                            function(x) all(x[!zv] == 0)), ,
                                     drop=FALSE]
        }
        hyp.matrix
    }

    if (ddf != "asymptotic") {
        ddf_info <- .Anova_ddf_setup(mod, component, ddf, not.aliased)
        rows <- .rbind_terms(
            lapply(seq_len(n.terms), function(term) {
                .Anova_ddf_test(mod, term.hyp.matrix(term), ddf, ddf_info)
            }),
            c("Fstat", "ndf", "ddf", "pval"))
        return(.Anova_ddf_result_table(names, ddf, rows, type = "III",
                                       response = responseName.default(mod)))
    }

    p <- teststat <- df <- rep(0, n.terms)
    for (term in seq_len(n.terms)){
        hyp.matrix <- term.hyp.matrix(term)
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
    }
    result <- data.frame(teststat, df, p)
    row.names(result) <- names
    names(result) <- c ("Chisq", "Df", "Pr(>Chisq)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c("Analysis of Deviance Table (Type III Wald chisquare tests)\n",
                                 paste("Response:", responseName.default(mod)))
    result
}
