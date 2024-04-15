
prior_synonyms <- c("fixef" = "beta",
                    "fixef_zi" = "beta_zi",
                    ## FIXME: update to betadisp when RE_disp is merged
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

    ## set up prior vectors (FIXME: should be a list?)
    for (p in prior_ivars) {
        assign(p, integer(np))
    }

    prior_params <- list()

    ## loop over prior specifications
    for (i in seq_len(np)) {

        ## process prior value (character to distribution code and parameter vector)
        pp <- priors[["prior"]][i]
        pname <- gsub("\\(.*", "", pp) ## extract distribution name
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
        ## if blank, all (except intercept)
        ## if names, locate the corresponding element(s)
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

        ## set pc to blank if missing
        pc <- priors[["coef"]][i] %||% "" |> trimws()
        if (pc == "") {
            if (substr(cl, 1, 4) == "beta") {
                prior_elend[i] <- length(info$fix[[match_names(cl)]])-1
                if ("(Intercept)" %in% info$fix[[match_names(cl)]]) {
                    prior_elstart[i] <- 1 ## skip intercept
                    ## (assume intercept is first column of model matrix ...)
                }
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
                    ## TO DO: match SD/cor components *within* log-SD vec??
                    component <- match_names(cl, prefix = "theta")
                    re_info <- info$re[[component]]
                    re_term <- nospace(pc)
                    ## matching entire component: try to match to full term (space-squashed)
                    w <- match(re_term, nospace(names(re_info)))
                    if (is.na(w)) {
                        ## ... or just grouping variable
                        w <- match(re_term, gsub("^[^|]+\\|", "", nospace(names(re_info))))
                    }
                    if (is.na(re_term)) stop("can't match prior RE term", re_term)
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
                            ## all suffixes
                            prior_elstart[i] <- theta_start
                            prior_elend[i] <- theta_start + nsd - 1
                        } else if (suffix == "cor") {
                            if (nsd == re_info[[w]]$blockNumTheta) {
                                warning("RE term has no correlation parameters")
                                ## FIXME: is it dangerous to proceed in this case?
                                ## stop/cancel out prior spec?
                            }
                            prior_elstart[i] <- theta_start + nsd
                            prior_elend[i] <- theta_start + re_info[[w]]$blockNumTheta - 1
                        } else stop("unknown prior suffix ", suffix)
                   } ## suffix specified 
                } ## specified theta elements
            } ## specified elements
        } ## component specified

        ## FIXME: don't hard-code here?
        prior_npar[i] <- switch(pname,
                                normal =,
                                gamma =,
                                cauchy =,
                                beta = 2,
                                t = 3,
                                lkj = 1,
                                other = stop("unknown prior type")
                                )

        if ((np <- length(prior_params[[i]])) != prior_npar[i]) {
            stop(sprintf("incorrect number of parameters for prior %s distribution (%d, was expecting %d)",
              pname, np, prior_npar[i]))
        }
        
    } ## loop over priors
    prior_params <- if (np == 0) numeric(0) else unlist(prior_params)
    ## FIXME: replace with mget() ?
    return(namedList(prior_distrib, prior_whichpar, prior_elstart, prior_elend, prior_npar, prior_params))
}

#' use of priors in glmmTMB
#'
#' (EXPERIMENTAL/subject to change)
#' 
#' \code{glmmTMB} can accept prior specifications, for doing maximum \emph{a posteriori} (MAP) estimation (or Hamiltonian MC with the \code{tmbstan} package), or (outside of a Bayesian framework) for the purposes of regularizing parameter estimates  
#' 
#' The \code{priors} argument to \code{glmmTMB} must (if not NULL) be a data frame with columns
#' \describe{
#' \item{\code{prior}}{character; the prior specification, e.g. "normal(0,2)"}
#' \item{\code{class}}{the name of the underlying parameter vector on which to impose the prior ("fixef", "fixef_zi", "fixef_disp", "ranef", "ranef_zi", "psi")}
#' \item{\code{coef}}{(optional) a string (if present) specifying the particular elements of the parameter vector to apply the prior to. \code{coef} should specify an integer parameter index, a column name from the fixed effect model matrix or a grouping variable for a random effect (the behaviour is currently undefined if there is more one than random effect term with the same grouping variable in a model ...); one can also append "_cor" or "_sd" to a random-effects \code{class} specification to denote the correlation parameters, or all of the standard deviation parameters, corresponding to a particular random effect term. If the \code{class} element is missing, or a particular element is blank, then all of the elements of the specified parameter vector use independent priors with the given specification. The exception is for fixed-effect parameter vectors, where the intercept (if present) is not included; the prior on the intercept must be set explicitly.}
#' }
#' `The available prior distributions are:
#' \itemize{
#' \item "normal" (mean/sd parameterization)
#' \item "t" (mean/sd/df)
#' \item "cauchy" (location/scale)
#' \item "gamma" (mean/shape); applied on the SD (\emph{not} the log-SD) scale
#' \item "lkj" (correlation) [WARNING, maybe buggy at present!]
#' }
#' The first three are typically used for fixed effect parameters; the fourth for standard deviation parameters; and the last for correlation structures. See the "priors" vignette for examples and further information.
#' 
#' @name priors
#'
#' @examples
#'
#' data("sleepstudy", package = "lme4")
#' prior1 <- data.frame(prior = c("normal(250,3)","t(0,3,3)","gamma(10,1)"),
#'                      class = c("fixef", "fixef", "ranef_sd"),
#'                      coef = c("(Intercept)", "Days", "Subject"))
#' g1 <- glmmTMB(Reaction ~ 1 + Days + (1 + Days |Subject), sleepstudy)
#' update(g1, prior = prior1)
#' prior2 <- data.frame(prior = c("t(0,3,3)","gamma(10,1)"),
#'                      class = c("fixef", "ranef_sd"),
#'                      coef = c("", "Subject"))
#' update(g1, prior = prior2) 
NULL


#' @importFrom stats reformulate
#' @importFrom utils capture.output
print.glmmTMB_prior <- function(x, compact = FALSE, ...) {
    if (is.null(x)) return(invisible(x))
    pstr <- character(nrow(x))
    for (i in seq_len(nrow(x))) {
        resp <- from_prior_syn(x$class[i])
        if (nzchar(x$coef[i])) {
            resp <- sprintf("%s (%s)", resp, x$coef[i])
        }
        ff <- reformulate(x$prior[i], response = resp)
        pstr[i] <- capture.output(print(showEnv = FALSE, ff))
        if (!compact) {
            cat(pstr[i], "\n")
        }
    }
    if (compact) cat(paste(pstr, collapse = "; "), "\n")
    invisible(x)
}
             
