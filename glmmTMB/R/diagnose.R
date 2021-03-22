##' diagnose model problems
##'
##' \strong{EXPERIMENTAL}. For a given model, this function attempts to isolate
##' potential causes of convergence problems. It checks (1) whether there are
##' any unusually large coefficients; (2) whether there are any unusually
##' scaled predictor variables; (3) if the Hessian (curvature of the
##' negative log-likelihood surface at the MLE) is positive definite
##' (i.e., whether the MLE really represents an optimum). For each
##' case it tries to isolate the particular parameters that are problematic.
##'
##' Problems in one category (e.g. complete separation) will generally
##' also appear in "downstream" categories (e.g. non-positive-definite
##' Hessians).
##' Therefore, it is generally
##' advisable to try to deal with problems in order, e.g. address problems with
##' complete separation first, then re-run the diagnostics to see whether
##' Hessian problems persist.
##' 
##' @param fit a \code{glmmTMB} fit
##' @param eval_eps numeric tolerance for 'bad' eigenvalues
##' @param evec_eps numeric tolerance for 'bad' eigenvector elements
##' @param big_coef numeric tolerance for large coefficients
##' @param big_sd_log10 numeric tolerance for badly scaled parameters (log10 scale), i.e. for default value of 3, predictor variables with sd less than 1e-3 or greater than 1e3 will be flagged)
##' @param big_zstat numeric tolerance for Z-statistic 
##' @param check_coefs identify large-magnitude coefficients? (Only checks conditional-model parameters if a (log, logit, cloglog, probit) link is used. Always checks zero-inflation, dispersion, and random-effects parameters. May produce false positives if predictor variables have extremely large scales.)
##' @param check_hessian identify non-positive-definite Hessian components?
##' @param check_zstats identify parameters with unusually large Z-statistics (ratio of standard error to mean)? Identifies likely failures of Wald confidence intervals/p-values.
##' @param check_scales identify predictors with unusually small or large scales?
##' @return a logical value based on whether anything questionable was found
##' @importFrom numDeriv jacobian hessian
##' @importFrom stats sd
##' @export
##'
diagnose <- function(fit, eval_eps=1e-5,evec_eps=1e-2,
                     big_coef=10,
                     big_sd_log10=3,
                     big_zstat=5,
                     check_coefs=TRUE,
                     check_zstats=TRUE,
                     check_hessian=TRUE,
                     check_scales=TRUE) {
    model_OK <- TRUE
    ## pull out the TMB object from the fit
    obj <- fit$obj
    ee <- environment(obj$fn)
    ## extract parameters
    pp <- ee$last.par[-ee$random]
    ss <- summary(fit$sdr)
    ss <- ss[grepl("^(beta|theta)",rownames(ss)),]
    ## easiest way to get names corresponding to all of the parameters
    nn <- tryCatch(colnames(vcov(fit,full=TRUE)),
                   ## fall-back position
                   error = function(e) make.unique(names(pp)))
    nn0 <- names(pp)
    names(pp) <- rownames(ss) <- nn
    ## check coefficients
    if (check_coefs) {
        link_par <- (nn0 != "beta" |
                     family(fit)$link %in% c("log","cloglog","logit","probit"))
        bigcoef <- (pp[abs(pp)>big_coef & link_par])
        if (length(bigcoef)>0) {
            model_OK <- FALSE
            cat(sprintf("Unusually large coefficients (|x|>%g):\n\n",big_coef))
            print(bigcoef)
            cat("\n")
            cat(strwrap(paste("Large negative coefficients in zi (log-odds of zero-inflation), dispersion, or random effects (log-standard deviations) suggest",
                              "unnecessary components (converging to zero on the constrained scale); large negative and/or positive",
                              "components in binomial or Poisson conditional parameters suggest (quasi-)complete separation",
                              collapse="")),"\n",sep="\n")
        }
    } ## check_coefs
    if (check_scales) {
        all_X <- lapply(c("X","Xzi","Xd"), function(x) getME(fit,x))
        all_X <- do.call(cbind,all_X)
        colnames(all_X) <- nn[seq(ncol(all_X))]
        sdvec <- apply(all_X, 2, sd)
        ## zero-variance columns excluded (presumably intercepts,
        ##  if not will presumably be caught elsewhere as collinear with
        ##  with intercepts.  Detecting "(Intercept)" seems too fragile?
        sdvec <- sdvec[sdvec>0 & abs(log10(sdvec))>big_sd_log10]
        if (length(sdvec)>0) {
            model_OK <- FALSE
            cat(sprintf("\npredictors with unusually large or small standard deviations (|log10(sd)|>%g):\n\n",sdvec))
            print(sdvec)
            cat("\n",strwrap(paste("Predictor variables with very narrow or wide ranges generally give rise to parameters with very large or",
                                   "small magnitudes, which can sometimes exacerbate numerical instability, and may also be appear",
                                   "(incorrectly) to be indicating a poorly defined optimum (i.e., a non-positive definite Hessian",
                                   collapse=" ")),sep="\n")
        }
    }
    if (check_zstats) {
        z <- ss[,"Std. Error"]/ss[,"Estimate"]
        bigz <- z[abs(z)>big_zstat]
        if (length(bigz)>0) {
            model_OK <- FALSE
            cat(sprintf("Unusually large Z-statistics (|x|>%g):\n\n",big_zstat))
            print(bigz)
            cat("\n")
            cat(strwrap(paste("Large Z-statistics (estimate/std err) suggest a failure ",
                              "of the Wald approximation - often also associated with ",
                              "parameters that are at or near the edge of their range ",
                              "(e.g. random-effects standard deviations approaching 0). ",
                              "While the Wald p-values and standard errors listed in ",
                              "summary() are unreliable, profile confidence intervals ",
                              "(see ?confint.glmmTMB) and likelihood ratio test p-values ",
                              "derived by comparing models (e.g. ?drop1) may still be OK. ",
                              "(Note that the LRT is conservative when the null value is ",
                              "on the boundary, e.g. a variance or zero-inflation value of 0 ",
                              "(Self and Liang 1987; Stram and Lee 1994; Goldman and Whelan 2000); ",
                              "in simple cases the p-value is approximately twice as large as it should be.)",
                              collapse="")),"\n",sep="\n")
        }
    }
    if (check_hessian) {
        if (!"sdr" %in% names(fit)) {
            warning("sdreport was not computed, skipping Hessian check")
        }
        if ("sdr" %in% names(fit) && !fit$sdr$pdHess) {
            model_OK <- FALSE
            cat("Non-positive definite Hessian\n\n")
            cat("parameters with non-finite standard deviations:\n")
            cat(strwrap(paste(nn[!is.finite(suppressWarnings(sqrt(diag(fit$sdr$cov.fixed))))],
                              collapse=", ")),"\n",sep="\n")
            h <- numDeriv::jacobian(obj$gr, pp)
            ## FIXME: consider SVD?
            eigs <- eigen(h)
            ## non-positive definite means some of the eigenvectors are <= 0
            bad <- which(eigs$values/max(eigs$values)<=eval_eps)
            if (length(bad)==0) {
                cat("Hessian seems OK\n")
                return(invisible(h))
            }
            cat(sprintf("maximum Hessian eigenvalue = %1.3g",eigs$values[1]),"\n")
            ## there could be more than one 'bad' direction/eigenvector ..
            for (b in bad) {
                cat(sprintf("Hessian eigenvalue %d = %1.3g (relative val = %1.3g)",
                            b,eigs$values[b],eigs$values[b]/eigs$values[1]),"\n")
                bad_vec <- eigs$vectors[,b]
                bad_elements <- which(abs(bad_vec)>evec_eps)
                cat("   bad elements:",nn[bad_elements],"\n")
            }
        } ## bad hessian
    } ## check hessian
    if (model_OK) cat("model looks OK!\n")
    return(invisible(model_OK))
}
