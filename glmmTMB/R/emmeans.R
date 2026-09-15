## methods for extending emmeans to handle glmmTMB objects

## NOTE: methods are dynamically exported by emmeans utility -- see code in zzz.R

##' Downstream methods
##'
##' @name downstream_methods
##' @aliases emmeans.glmmTMB
##'
##' @description
##' Methods have been written that allow \code{glmmTMB} objects to be used with
##' several downstream packages that enable different forms of inference.
##' For some methods (\code{Anova} and \code{emmeans}, but \emph{not} \code{effects} at present),
##' set the \code{component} argument
##' to "cond" (conditional, the default), "zi" (zero-inflation) or "disp" (dispersion) in order to produce results
##' for the corresponding part of a \code{glmmTMB} model. 
##' Support for \pkg{emmeans} also allows additional options 
##' \code{component = "response"} (response means taking both the \code{cond} and
##' \code{zi} components into account), and \code{component = "cmean"} (mean of the 
##' [possibly truncated] conditional distribution). 
##'
##' In particular,
##' \itemize{
##' \item \code{car::Anova} constructs type-II and type-III Anova tables
##' for the fixed effect parameters of any component
##' \item the \code{emmeans} package computes estimated marginal means (previously known as least-squares means)
##' for the fixed effects of any component, or predictions with \code{type = "response"} or
##' \code{type = "component"}. Note: In hurdle models, 
##' \code{component = "cmean"} produces means
##' of the truncated conditional distribution, while 
##' \code{component = "cond", type = "response"} produces means of the \emph{untruncated}
##' conditional distribution.
##' \item the \code{effects} package computes graphical tabular effect displays
##' (only for the fixed effects of the conditional component)
##' }
##' @section Denominator degrees of freedom in \code{emmeans}:
##' For models with random effects, the \code{ddf} argument to \code{emmeans()}
##' (default taken from \code{getOption("glmmTMB.df", "asymptotic")}) additionally accepts
##' \code{"satterthwaite"} and \code{"kenward-roger"} (see \code{\link{dof_KR}} and
##' \code{\link{dof_satt}} for the underlying calculations), matching the same argument
##' to \code{\link{summary.glmmTMB}} and \code{\link{anova.glmmTMB}}. \code{ddf = "kenward-roger"}
##' requires a model fitted with \code{REML = TRUE}: for an ML fit (\pkg{glmmTMB}'s default),
##' it throws an error rather than silently substituting another method; it also requires a
##' family with an estimated dispersion parameter, and throws an error for families such as
##' \code{binomial} or \code{poisson} that lack one. For families other than \code{gaussian},
##' \code{"kenward-roger"} and \code{"satterthwaite"} are allowed but emit a warning, because
##' their performance (and theoretical justification) for GLMMs is poorly understood.
##'
##' For Gaussian models \emph{without} random effects, \code{emmeans()} defaults
##' to the residual degrees of freedom for a plain fit, i.e. one with
##' \code{dispformula = ~1} and an estimated dispersion parameter. Any other such
##' fit defaults to \code{"asymptotic"} (infinite df): a non-trivial
##' \code{dispformula}, \code{dispformula = ~0}, or a dispersion parameter held
##' fixed via the \code{map} argument to \code{\link{glmmTMB}}. In the last case
##' there is no variance parameter left to estimate, so the residual variance is
##' known and the Wald statistics are exactly standard normal. (As elsewhere in
##' \pkg{glmmTMB}, the residual degrees of freedom count the dispersion
##' parameter, so they are one lower than \code{lm()} reports for the same
##' fixed-effect model.) Those are defaults, as is a value taken from
##' \code{getOption("glmmTMB.df")}; a \code{ddf} passed in the call itself is
##' respected where possible. \code{"kenward-roger"} and \code{"satterthwaite"}
##' need random effects, so for models without them they fall back to the
##' residual degrees of freedom with a message, exactly as in
##' \code{\link{summary.glmmTMB}}. That fallback also applies to a model whose
##' dispersion parameter is fixed, where it overrides the infinite-df default
##' described above, so that \code{emmeans()} and \code{summary()} give the same
##' answer to the same request.
##' @param mod a glmmTMB model
##' @param component which component of the model to test/analyze ("cond", "zi", or "disp")
##'     or, in \pkg{emmeans} only, "response" or "cmean" as described in Details.
##' @param \dots Additional parameters that may be supported by the method.
##' @details While the examples below are disabled for earlier versions of
##' R, they may still work; it may be necessary to refer to private
##' versions of methods, e.g. \code{glmmTMB:::Anova.glmmTMB(model, ...)}.
##' @importFrom stats delete.response
##' @examples
##' warp.lm <- glmmTMB(breaks ~ wool * tension, data = warpbreaks)
##' salamander1 <- up2date(readRDS(system.file("example_files","salamander1.rds",package="glmmTMB")))
##' if (require(emmeans)) withAutoprint({
##'     emmeans(warp.lm, poly ~ tension | wool)
##'     emmeans(salamander1, ~ mined, type="response")  # conditional means
##'     emmeans(salamander1, ~ mined, component="cmean")     # same as above, but re-gridded
##'     emmeans(salamander1, ~ mined, component="zi", type="response")  # zero probabilities
##'     emmeans(salamander1, ~ mined, component="response")  # response means including both components
##' })
##' if (getRversion() >= "3.6.0") {
##'    if (require(car)) withAutoprint({
##'        Anova(warp.lm,type="III")
##'        Anova(salamander1)
##'        Anova(salamander1, component="zi")
##'    })
##'    if (require(effects)) withAutoprint({
##'        plot(allEffects(warp.lm))
##'        plot(allEffects(salamander1))
##'    })
##' }
NULL  ## don't document the files here!


## recover_data method -- DO NOT export -- see zzz.R
## do not document either

recover_data.glmmTMB <- function (object, component = c("cond", "zi", "disp", "response", "cmean"), ...) {
    if (!requireNamespace("emmeans")) 
        stop("please install (if necessary) and load the emmeans package")
    
    component <- match.arg(component)
    # which terms to use?
    tcomp <- ifelse(component %in% c("response", "cmean"), "cond", component)
    trms <- delete.response(terms(object, component = tcomp))
    nbasis <- estimability::all.estble
    if (component %in% c("response", "cmean")) {  # may need add'l terms for response mode
        if(!is.null(ztrms <- terms(object, component = "zi")) && (length(all.vars(ztrms)) > 0))
            trms <- emmeans::.combine.terms(trms, ztrms)
        if(!is.null(dtrms <- terms(object, component = "disp")) && (length(all.vars(dtrms)) > 0))
            trms <- emmeans::.combine.terms(trms, dtrms)
    }
    else if (component != "cond") {
        if (is.null(trms) || (length(all.vars(trms)) == 0))
            stop("No reference grid is available for the '", component, "' component")
    }
    fcall <- getCall(object)
    emmeans::recover_data(fcall, trms, 
                          attr(model.frame(object), "na.action"), ...)
}


emm_basis.glmmTMB <- function (object, trms, xlev, grid, component = c("cond", "zi", 
                                                                       "disp", "response", "cmean"), vcov.,
                               ddf = getOption("glmmTMB.df", "asymptotic"),  ...) {

    ## FIXME: implement a 'KR limit' argument/option that determines whether to use KR for large problems ... ??
    component <- match.arg(component)
    check_dots(.ignore = c("misc", "options"))
    misc <- list()
    ## ddf-processing
    ## 1. no random effects
    fam <- family(object)$family

    ## did the caller actually ask for a particular ddf in this call, or are we
    ## falling back on the default? only the latter may be silently overridden
    ## below. A ddf coming from getOption("glmmTMB.df") is a default, not a
    ## request: summary()/anova()/Anova() don't read that option at all, so
    ## treating it as a request here would make emmeans() disagree with them
    ## whenever it is set. NB has to be evaluated before the match.arg() below,
    ## which would make missing(ddf) FALSE
    ddf_explicit <- !missing(ddf)
    ## same choices as summary()/anova()/Anova(); without this an unrecognized
    ## string silently ended up as residual df. "df.residual" is what get_ddf()
    ## returns internally and has always been accepted here too, so keep it
    ## working, but don't offer it as a choice: the other entry points don't
    ## take it
    if (!identical(ddf, "df.residual")) {
        ddf <- match.arg(ddf, c("asymptotic", "kenward-roger", "satterthwaite"))
    }

    ddf_set <- function(used, requested = ddf) {
        if (requested != used) {
            warning(gettextf("ddf '%s' specified, using ddf '%s' instead", requested, used))
        }
        return(used)
    }
    get_ddf <- function() {

        if (component != "cond") return(ddf_set("asymptotic"))

        if (!hasRandom(object)) {
            if (fam != "gaussian") return(ddf_set("asymptotic"))
            if (!ddf_explicit) {
                ## default: residual df for a plain LM-like fit (nobs - npar,
                ## which counts the dispersion parameter, so one fewer than
                ## lm() would report) -- but *not* when no dispersion
                ## parameter is estimated at all (e.g. pinned via 'map'). The
                ## residual variance is known then, the Wald statistics are
                ## exactly normal, and residual df would only make the
                ## intervals spuriously wide
                if (trivialDisp(object) && estDisp(object)) {
                    return("df.residual")  ## don't want to warn here
                }
                return("asymptotic")
            }
            ## an explicit request is honoured where it can be, and otherwise
            ## downgraded through the same check_ddf() that
            ## summary()/anova()/Anova() use, so emmeans gives the same
            ## message and the same df instead of silently ignoring the
            ## request (trivial dispformula) or erroring inside
            ## GMRFmarginal() on its way through dof_satt() (non-trivial
            ## one), which needs the joint precision matrix of a model with
            ## random effects
            if (ddf %in% c("asymptotic", "df.residual")) return(ddf)
            check_ddf(object, ddf)
            return("df.residual")
        }

        ## hard error (not a silent downgrade) for KR + non-REML, matching
        ## summary.glmmTMB()/anova.glmmTMB() via check_ddf()
        if (ddf == "kenward-roger") {
            .check_KR_reml(object)
            ## Kenward-Roger's variance-component machinery only supports
            ## families with an estimated dispersion parameter (see the
            ## matching check in check_ddf()); without this, families such
            ## as binomial/poisson fail with an opaque error instead
            if (!usesDispersion(fam)) {
                stop(sprintf(
                    "ddf='kenward-roger' is not supported for family '%s' (no estimated dispersion parameter); use ddf='satterthwaite' or ddf='asymptotic' instead",
                    fam), call. = FALSE)
            }
        }

        if (fam != "gaussian" && ddf != "asymptotic") .warn_ddf_glmm(ddf)
        return(ddf)
    }

    ddf <- get_ddf()

    if (ddf == "kenward-roger") {
        V <- vcov(object)[[component]]
        dfargs <- list(unadjV = V,
                       adjV = .vcov_kenward_adjusted(object))
        V_kr <- as.matrix(dfargs$adjV)
        V <- V_kr
        dffun <- function(k, dfargs) pbkrtest::Lb_ddf(k, dfargs$unadjV, dfargs$adjV)
    } else if (ddf == "satterthwaite") {
        ## emmeans::ref_grid() strips dffun's enclosing environment
        ## (sets it to baseenv()), so dffun can't rely on free variables
        ## such as a captured dof_satt (glmmTMB#1304) -- stash it in
        ## dfargs instead, where it's reached via the 'dfargs' argument
        dfargs <- list(object = object, dof_satt = dof_satt)
        ## emmeans calls dffun() once per contrast, passing a bare
        ## vector k rather than a full contrast matrix; dof_satt()
        ## expects a matrix (one row per contrast), so wrap k accordingly
        dffun <- function(k, dfargs) suppressMessages(dfargs$dof_satt(dfargs$object, L = matrix(k, nrow = 1)))
    } else if (ddf == "df.residual") {
        dfargs <- list(object = object)
        dffun <- function(k, dfargs) stats::df.residual(dfargs$object)
    } else if (ddf == "asymptotic") {
        dfargs <- list()
        dffun <- function(k, dfargs) Inf
    } else stop(sprintf("unknown ddf specification '%s'", ddf))

    # internal fcn for identifying non-estimable components
    .which.nonest <- function(cmp) {
        bh <- fixef(object)[[cmp]]
        if (!any(is.na(bh)))
            return(numeric(0))  # no estimability issues
        tms <- delete.response(terms(object, component = cmp))
        bas <- emm_basis.glmmTMB(object, tms, xlev, grid, component = cmp)
        which(!estimability::is.estble(bas$X, bas$nbasis))
    }
    
    nbasis <- estimability::all.estble
    if (component %in% c("response", "cmean")) {
        ptype <- ifelse(component == "cmean", "conditional", 
                        "response")
        for (nm in object$modelInfo$grpVar) grid[[nm]] <- NA
        tmp <- predict(object, newdata = grid, type = ptype, 
                       re.form = NA, se.fit = TRUE, cov.fit = TRUE)
        bhat <- tmp$fit
        X <- diag(1, length(bhat))
        V <- tmp$cov.fit
        if(component == "response")
            bhat[.which.nonest("zi")] <- NA
        bhat[.which.nonest("cond")] <- NA
        bhat[.which.nonest("disp")] <- NA
        if (length(w <- which(is.na(bhat))) > 0) {
            nbasis <- matrix(0, nrow = length(bhat), ncol = length(w))
            for (j in seq_along(w))
                nbasis[w[j], j] <- 1
            V <- V[-w, -w, drop = FALSE]
        }
    }
    else {
	## combinomial with allow_negative_nu uses identity link on dispersion;
        ## other families always use log link on dispformula
        disp_link <- if (isTRUE(object$modelInfo$family$allow_negative_nu)) "identity" else "log"
        fam <- switch(component, cond = family(object), zi = list(link = "logit"), 
                      disp = list(link = disp_link))
        misc <- emmeans::.std.link.labels(fam, misc)
        if (missing(vcov.)) {
            V <- as.matrix(vcov(object, include_nonest = FALSE)[[component]])
            ## coefficients fixed via 'map' are known constants: pad the
            ## covariance matrix with zero rows/columns so its dimension
            ## matches the full coefficient vector used for the grid
            V <- pad_mapped_vcov(object, V, component)
        }
        else {
            V <- vcov.
        }
        contrasts <- attr(model.matrix(object, component = component), 
                          "contrasts")
        m <- model.frame(trms, grid, na.action = na.pass, xlev = xlev)
        X <- model.matrix(trms, m, contrasts.arg = contrasts)
        bhat <- fixef(object)[[component]]
        if(any(is.na(bhat))) {
            modmat <- model.matrix(trms, model.frame(object), 
                                   contrasts.arg = contrasts)
            nbasis <- estimability::nonest.basis(modmat)
        }
    }
    namedList(X, bhat, nbasis, V, dffun, dfargs, misc)
}

