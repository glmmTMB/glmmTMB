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
##' @param explain provide detailed explanation of each test?
##' @return a logical value based on whether anything questionable was found
##' @importFrom numDeriv jacobian hessian
##' @importFrom stats sd
##' @export
##'
diagnose <- function(fit,
                     eval_eps=1e-5,
                     evec_eps=1e-2,
                     big_coef=10,
                     big_sd_log10=3,
                     big_zstat=5,
                     check_coefs=TRUE,
                     check_zstats=TRUE,
                     check_hessian=TRUE,
                     check_scales=TRUE,
                     explain = TRUE) {
    prt_explain <- function(...) {
      if (explain) {
        s <- do.call(paste, c(list(...), list(collapse = "")))
        cat(strwrap(s), "\n", sep = "\n")
      }
      return(invisible(NULL))
    }
    model_OK <- TRUE
    ## pull out the TMB object from the fit
    obj <- fit$obj
    ee <- obj$env
    ## extract parameters
    pp <- ee$last.par.best
    if (!is.null(r <- ee$random)) { pp <- pp[-r] }
    ss <- suppressWarnings(summary(fit$sdr))
    ss <- ss[grepl("^(beta|theta)", rownames(ss)), ]
    ## easiest way to get names corresponding to all of the parameters
    vv <- try(vcov(fit, full = TRUE))
    if (!inherits(vv, "try-error")) {
        nn <- colnames(vcov(fit, full = TRUE))
    } else {
        nn <- make.unique(names(pp))
        if (isREML(fit)) warning("bad vcov/REML combination is not tested, may fail")
    }
    nn0 <- names(pp)
    if (isREML(fit)) {
        ss <- ss[rownames(ss) != "beta",]
        nn <- nn[!startsWith(names(nn), "cond")]
    }
    rownames(ss) <- names(pp) <- nn
    ## check coefficients
    if (check_coefs) {
        ## logic: if link function for *conditional* model is unitless then all parameters need checking
        ## otherwise only check parameters *other* than conditional model parameters ("betad", "betazi", "theta")
        link_par <- (nn0 != "beta" |
                     family(fit)$link %in% c("log", "cloglog", "logit", "probit"))
        bigcoef <- (pp[abs(pp)>big_coef & link_par])
        if (length(bigcoef)>0) {
            model_OK <- FALSE
            cat(sprintf("Unusually large coefficients (|x|>%g):\n\n",big_coef))
            print(bigcoef)
            cat("\n")
            ## FIXME: could do more for the user by checking
            ## for the user *which* elements were large/small
            ## rather than giving them the explanation of all possible cases
            prt_explain("Large negative coefficients in zi (log-odds of zero-inflation), dispersion, or random effects (log-standard deviations) suggest",
                        "unnecessary components (converging to zero on the constrained scale); large negative and/or positive",
                        "components in binomial or Poisson conditional parameters suggest (quasi-)complete separation.",
                        "Large values of nbinom2 dispersion suggest that you should use a Poisson model instead.")
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
            prt_explain("Predictor variables with very narrow or wide ranges generally give rise to parameters with very large or",
                        "small magnitudes, which can sometimes exacerbate numerical instability, and may also be appear",
                        "(incorrectly) to be indicating a poorly defined optimum (i.e., a non-positive definite Hessian")
        }
    }
    if (check_zstats) {
        z <- ss[,"Estimate"]/ss[,"Std. Error"]
        bigz <- z[!is.na(z) & abs(z)>big_zstat]
        if (length(bigz)>0) {
            model_OK <- FALSE
            cat(sprintf("Unusually large Z-statistics (|x|>%g):\n\n",big_zstat))
            print(bigz)
            cat("\n")
            prt_explain("Large Z-statistics (estimate/std err) suggest a *possible* failure ",
                        "of the Wald approximation - often also associated with ",
                        "parameters that are at or near the edge of their range ",
                        "(e.g. random-effects standard deviations approaching 0).",
                        " (Alternately, they may simply represent very well-estimated parameters;",
                        "intercepts of non-centered models may fall in this category.)",
                        "While the Wald p-values and standard errors listed in ",
                        "summary() may be unreliable, profile confidence intervals ",
                        "(see ?confint.glmmTMB) and likelihood ratio test p-values ",
                        "derived by comparing models (e.g. ?drop1) are probably still OK. ",
                        "(Note that the LRT is conservative when the null value is ",
                        "on the boundary, e.g. a variance or zero-inflation value of 0 ",
                        "(Self and Liang 1987; Stram and Lee 1994; Goldman and Whelan 2000); ",
                        "in simple cases these p-values are approximately twice as large as they should be.)")
        }
    }
    if (check_hessian) {
        if (!"sdr" %in% names(fit)) {
            warning("sdreport was not computed, skipping Hessian check")
        }
        if ("sdr" %in% names(fit) && !fit$sdr$pdHess) {
            model_OK <- FALSE
            cat("Non-positive definite (NPD) Hessian\n\n")
            prt_explain("The Hessian matrix represents the curvature of the",
                        "log-likelihood surface at the maximum likelihood estimate (MLE) of the parameters",
                        "(its inverse is the estimate of the parameter covariance matrix). ",
                        "A non-positive-definite Hessian means that the likelihood surface is approximately flat ",
                        "(or upward-curving) at the MLE, which means the model is overfitted or poorly posed ",
                        "in some way.",
                        "NPD Hessians are often associated with extreme parameter estimates."

                        )
            nonfinite_sd <- !is.finite(suppressWarnings(sqrt(diag(fit$sdr$cov.fixed))))
            if (any(nonfinite_sd)) {
              cat("parameters with non-finite standard deviations:\n")
              if (all(nonfinite_sd)) {
                cat("(all of them!)\n\n")
              } else {
                cat(strwrap(paste(nn[nonfinite_sd],
                                  collapse=", ")),"\n\n",sep="\n")
              }
            }
            ## fit hessian with Richardson extrapolation (more accurate/slower than built-in optimHess)
            cat("recomputing Hessian via Richardson extrapolation.",
                "If this is too slow, consider setting check_hessian = FALSE",
                "\n\n")
            h <- numDeriv::jacobian(obj$gr, pp)
            ## FIXME: consider SVD?
            ## FIXME: add explanation
            eigs <- eigen(h)
            ## non-positive definite means some of the eigenvectors are <= 0
            complex_eigs <- is.complex(eigs$values)
            if (!complex_eigs) {
              bad <- which(eigs$values/max(eigs$values) <= eval_eps)
            }
            if (!complex_eigs && length(bad) == 0) {
              cat("Hessian seems OK\n")
              prt_explain("glmmTMB's internal calculations suggested that the Hessian was bad/non-positive definite;",
                          "however, a slower and more precise calculation suggests that it's actually OK. Your model",
                          "may be somewhat numerically unstable.")
              return(invisible(h)) ## bail out here
            }
            if (complex_eigs) {
              cat("Hessian has complex eigenvalues\n\n")
              prt_explain("We would have used the smallest eigenvalues of the ",
                          "Hessian to determine which components were bad",
                          "but instead we got complex eigenvalues.",
                          "(Not really sure what to do with this ...)")
            } else {
              prt_explain("The next set of diagnostics attempts to determine which elements of the Hessian",
                          "are causing the non-positive-definiteness. ",
                          "Components with very small eigenvalues represent 'flat' directions, ",
                          "i.e., combinations of parameters for which the data may contain very little information. ",
                          sprintf("So-called 'bad elements' represent the dominant components (absolute values >%1.3g)", evec_eps),
                          "of the eigenvectors corresponding to the 'flat' directions")
              cat(sprintf("maximum Hessian eigenvalue = %1.3g",eigs$values[1]),"\n")
              ## there could be more than one 'bad' direction/eigenvector ..
              for (b in bad) {
                cat(sprintf("Hessian eigenvalue %d = %1.3g (relative val = %1.3g)",
                            b,eigs$values[b],eigs$values[b]/eigs$values[1]),"\n")
                bad_vec <- eigs$vectors[,b]
                bad_elements <- which(abs(bad_vec)>evec_eps)
                cat("   bad elements:", nn[bad_elements], "\n")
            }
            } ## hessian has real eigs
        } ## bad hessian
    } ## check hessian
    if (model_OK) cat("model looks OK!\n")
    return(invisible(model_OK))
}
