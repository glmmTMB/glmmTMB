
prior_synonyms <- c("fixef" = "beta",
                    "fixef_zi" = "beta_zi",
                    "fixef_disp" = "betad",
                    "ranef" = "theta",
                    "ranef_zi" = "theta_zi",
                    "psi" = "shape")

to_prior_syn <- function(x) {
    if (!x %in% names(prior_synonyms)) return (x)
    prior_synonyms[x]
}

from_prior_syn <- function(x) {
    if (!x %in% prior_synonyms) return(x)
    names(prior_synonyms)[match(x, prior_synonyms)]
}
    
#' @noRd
#' @examples
#' if (require(brms)) {
#' bprior <- c(prior_string("normal(0,10)", class = "beta"),
#' ##            prior(normal(1,2), class = b, coef = treat),
#'            prior_(~cauchy(0,2), class = ~betad))
#' proc_priors(bprior)
#' }
proc_priors <- function(priors) {
    ## priors is a data frame as in brms
    ## process prior list into TMB data structures
    ##
    np <- if (is.null(priors)) 0 else nrow(priors)
    prior_distrib <- prior_whichpar <- prior_element <- integer(np)
    prior_params <- list()
    for (i in seq_len(np)) {
        ## prior is a string
        pp <- priors[["prior"]][i]
        pname <- gsub("\\(.*","",pp)
        prior_distrib[i] <- .valid_prior[pname]
        if (is.na(prior_distrib[i])) stop("unknown prior distribution ",pname)
        ## extract parameter values
        ## parse and drop expression()
        p_params <- parse(text = pp)[[1]]
        p_params[[1]] <- as.name("c")
        ## in which environment???
        prior_params[[i]] <- eval(p_params)
        cl <- to_prior_syn(priors[["class"]][i])
        prior_whichpar[i] <- .valid_vprior[cl]
        if (is.na(prior_whichpar[i])) stop("unknown prior variable ", cl)
                                           
        pc <- priors[["coef"]][i]
        if (pc == "") prior_element[i] <- NA_integer_
        if (pc != "") stop("element-specific priors not implemented yet")
    }
    prior_params <- if (np == 0) numeric(0) else unlist(prior_params)
    return(namedList(prior_distrib, prior_whichpar, prior_element,
                 prior_params))
}

#' use of priors in glmmTMB
#'
#' \code{glmmTMB} can accept prior specifications, for doing maximum a posteriori (MAP) estimation (or Hamiltonian MC with the \code{tmbstan} package), or outside of a strictly Bayesian framework for the purposes of regularizing estimates
#' 
#' The \code{priors} argument to \code{glmmTMB} must (if not NULL) be a data frame with columns (at least) \code{prior} (character; the prior specification, e.g. "normal(0,2)"); \code{class} (the name of the underlying parameter vector on which to impose the prior ("beta", "betazi", "betad", "theta", "thetazi", "psi"); \code{coef} a string specifying the particular elements of the parameter vector to apply the prior to (not yet implemented; must be specified as an empty string). At present priors can only be imposed jointly on all of the elements of a specified parameter vector, e.g. all fixed-effect coefficients. The tools in \code{brms} for specifying priors (e.g. \code{set_prior} should work to produce legal specifications. The available prior distributions are "normal" (mean/sd parameterization); "t" (mean/sd/df); "cauchy" (location/scale); "gamma" (mean/shape). The first three are typically used for fixed effect parameters; the last is typically used for standard deviation parameters ...
#'
#' to be continued ... (complete separation, singularity, etc.)
#' 
#' @name priors
#' 
NULL

