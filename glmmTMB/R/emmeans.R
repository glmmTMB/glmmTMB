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
                                                                       "disp", "response", "cmean"), vcov., ...) 
{
    component <- match.arg(component)
    L <- list(...)
    if (length(L) > 0) {
    }
    misc <- list()
    if (family(object)$family == "gaussian") {
        dfargs <- list(df = df.residual(object))
        dffun <- function(k, dfargs) dfargs$df
    }
    else {
        dffun <- function(k, dfargs) Inf
        dfargs <- list()
    }
    
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
        fam <- switch(component, cond = family(object), zi = list(link = "logit"), 
                      disp = list(link = "log"))
        misc <- emmeans::.std.link.labels(fam, misc)
        if (missing(vcov.)) {
            V <- as.matrix(vcov(object, include_nonest = FALSE)[[component]])
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

