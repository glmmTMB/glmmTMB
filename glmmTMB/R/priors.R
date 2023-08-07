
prior_synonyms <- c("fixef" = "beta",
                    "fixef_zi" = "beta_zi",
                    "fixef_disp" = "betad",
                    "ranef" = "theta",
                    "ranef_zi" = "theta_zi",
                    "psi" = "shape")

prior_ivars <- paste0("prior_", c("distrib", "whichpar", "elstart", "elend", "npar"))
prior_fvars <- "prior_params"

to_prior_syn <- function(x) {
    if (!x %in% names(prior_synonyms)) return (x)
    prior_synonyms[x]
}

from_prior_syn <- function(x) {
    if (!x %in% prior_synonyms) return(x)
    names(prior_synonyms)[match(x, prior_synonyms)]
}
    
#' @noRd
#' @param priors data frame with ...
#' @param info additional info needed: list(fix = list(cond = , zi =, disp = )) [colnames of X matrices];
#'  re = list(cond = list(cnms, ss), zi = ) [cnms, names of grouping vars and name of terms; ss, structure code]
#' @return a list with columns
#' \itemize{
#' \item prior_distrib: distribution code
#' \item prior_np: number of parameters
#' \item prior_whichpar: which parameter vector (theta, beta, etc.)
#' \item prior_elstart: starting element of parameter vector
#' \item prior_elend: ending element of parameter vector
#' \item prior_npars: number of prior (hyper)parameters for each term
#' \item prior_params: vector of all hyperparameters, concatenated
#' }
proc_priors <- function(priors, info = NULL) {
    ## priors is a data frame as in brms
    ## process prior list into TMB data structures
    ## 'info' parameter for translating elements into indices; not implemented yet
    np <- if (is.null(priors)) 0 else nrow(priors)
    ## FIXME: do this automatically via prior_ivars?
    prior_distrib <- prior_whichpar <- prior_elstart <- prior_elend <- prior_npar <- integer(np)
    prior_params <- list()
    for (i in seq_len(np)) {

        ## process prior value (character to distribution code and parameter vector)
        pp <- priors[["prior"]][i]
        pname <- gsub("\\(.*","",pp) ## strip prior name
        prior_distrib[i] <- .valid_prior[pname]
        if (is.na(prior_distrib[i])) stop("unknown prior distribution ",pname)
        
        ## extract parameter values
        ## parse and drop expression()
        p_params <- parse(text = pp)[[1]]
        p_params[[1]] <- as.name("c")
        ## in which environment???
        prior_params[[i]] <- eval(p_params)

        ## process 'class' (parameter vector)
        pcl <- priors[["class"]][i]
        if (!grepl("_(cor|sd)", pcl)) {
            suffix <- NA_character_
        } else {
            suffix <- gsub("^.*((cor|sd))$", "\\1", pcl)
        }
        pcl <- gsub("_(cor|sd)$", "", pcl)

        cl <- to_prior_syn(pcl)
        prior_whichpar[i] <- .valid_vprior[cl]
        if (is.na(prior_whichpar[i])) stop("unknown prior variable ", cl)

        ## figure out elements
        ## if blank, all
        ## if names, locate
        ## if non-blank suffix (sd/cor), figure out which elements based on ss/cnms
        ## process 'coef' (particular element)

        nthetavec <- sapply(info$re,
                            function(x) {
                                ntheta <- vapply(x, "[[",
                                                 "blockNumTheta",
                                                 FUN.VALUE = numeric(1))
                                cc <- cumsum(ntheta)
                                ## want *starting* value of each theta term
                                ## keep names, shift back one
                                nm <- names(cc)
                                if (length(cc) == 0) return(integer(0))
                                cc <- c(1, cc) |> setNames(c(nm, "..total"))
                                return(cc)
                            })

        nospace <- function(x) gsub(" +", "", x)
        thetanames <- lapply(info$re,
                             function(x) nospace(names(x)))
        ## 
        pc <- trimws(priors[["coef"]][i])
        if (pc == "") {
            if (substr(cl, 1, 4) == "beta") {
                prior_elend[i] <- length(info$fix[[match_names(cl)]])-1
            } else {
                prior_elend[i] <- nthetavec[[match_names(cl, prefix = "theta")]][["..total"]] -1 
            }
        } else {
            ## single numeric index (subtract 1 for R to C++ indexing shift)
            if (grepl("^[0-9]+$", pc)) {
                ## FIXME: should check vector length here (it is
                ## also caught inside TMB code)
                prior_elstart[i] <- prior_elend[i] <- as.integer(pc) -1
            } else {
                if (substr(cl, 1, 4) == "beta") {
                    ind <- match(pc, info$fix[[match_names(cl)]])
                    if (is.na(ind)) stop("can't match prior element ",pc)
                    prior_elstart[i] <- prior_elend[i] <- ind-1
                } else {
                    ## match component based on cnms
                    ## work out number of sd/cor params based on structure
                    ## first need to locate theta component in overall
                    ##  theta vector
                    component <- match_names(cl, prefix = "theta")
                    re_info <- info$re[[component]]
                    w <- match(nospace(pc), nospace(names(re_info)))
                    if (is.na(w)) stop("can't match prior RE component ", pc)
                    theta_start <- nthetavec[[component]][w] - 1 ## C++ index
                    if (is.na(suffix)) {
                        prior_elstart[i] <- theta_start
                        prior_elend[i] <- theta_start + re_info[[w]]$blockNumTheta - 1
                    } else {
                        blocksize <- re_info[[w]]$blockSize
                        blockcodelab <- names(.valid_covstruct)[match(re_info[[w]]$blockCode, .valid_covstruct)]
                        if (blockcodelab == "rr") stop("can't do priors for rr models yet")
                        nsd <- if (blockcodelab == "homdiag") 1 else blocksize
                        if (suffix == "sd") {
                            prior_elstart[i] <- theta_start
                            prior_elend[i] <- theta_start + nsd - 1
                        } else {
                            if (nsd == re_info[[w]]$blockNumTheta) {
                                warning("RE term has no correlation parameters")
                                ## FIXME: is it dangerous to proceed in this case?
                                ## stop/cancel out prior spec?
                            }
                            prior_elstart[i] <- theta_start + nsd
                            prior_elend[i] <- theta_start + re_info[[w]]$blockNumTheta - 1
                        }
                    }
                } ## specified theta elements
            } ## specified elements
        } ## component specified
        
        prior_npar[i] <- switch(pname,
                                normal =,
                                gamma =,
                                cauchy =,
                                beta = 2,
                                t = 3,
                                lkj = 1,
                                other = stop("unknown prior type")
                                )
        
    } ## loop over priors
    prior_params <- if (np == 0) numeric(0) else unlist(prior_params)
    ## FIXME: replace with mget() ?
    return(namedList(prior_distrib, prior_whichpar, prior_elstart, prior_elend, prior_npar, prior_params))
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

