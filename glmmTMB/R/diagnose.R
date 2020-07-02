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
##' @param big_sd numeric tolerance for badly scaled
##' @param check_coefs identify large-magnitude coefficients?
##' @param check_hessian identify non-positive-definite Hessian components?
##' @param check_scales identify predictors with unusually small or large scales?
##' @return a logical value based on whether anything questionable was found
##' @importFrom numDeriv jacobian hessian

##' @export
##'
diagnose <- function(fit, eval_eps=1e-5,evec_eps=1e-2,big_coef=10,big_log10_sd=3,
                     check_coefs=TRUE,
                     check_hessian=TRUE,
                     check_scales=TRUE) {
    model_OK <- TRUE
    ## pull out the TMB object from the fit
    obj <- fit$obj
    ee <- environment(obj$fn)
    ## extract parameters
    pp <- ee$last.par[-ee$random]
    ## easiest way to get names corresponding to all of the parameters
    nn <- tryCatch(colnames(vcov(fit,full=TRUE)),
                   ## fall-back position
                   error = function(e) make.unique(names(pp)))
    names(pp) <- nn
    ## check coefficients
    if (check_coefs) {
        if (length(bigcoef <- pp[abs(pp)>big_coef])>0) {
            model_OK <- FALSE
            cat(sprintf("Unusually large coefficients (|beta|>%g):\n\n",big_coef))
            print(bigcoef)
            cat("\n")
            cat(strwrap(paste("Large negative coefficients in zi (log-odds of zero-inflation) or random effects (log-standard deviations) suggest",
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
        sdvec <- sdvec[sdvec>0 & abs(log10(sdvec))>big_log10_sd]
        if (length(sdvec)>0) {
            model_OK <- FALSE
            cat(sprintf("\npredictors with unusually large or small standard deviations (|log10(sd)|>%g):\n\n",big_log10_sd))
            print(sdvec)
            cat("\n",strwrap(paste("Predictor variables with very narrow or wide ranges generally give rise to parameters with very large or",
                                   "small magnitudes, which can sometimes exacerbate numerical instability, and may also be appear",
                                   "(incorrectly) to be indicating a poorly defined optimum (i.e., a non-positive definite Hessian",
                                   collapse=" ")),sep="\n")
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
