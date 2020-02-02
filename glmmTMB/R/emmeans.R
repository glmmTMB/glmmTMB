## methods for extending emmeans to handle glmmTMB objects

## NOTE: methods are dynamically exported by emmeans utility -- see code in zzz.R

##' Downstream methods for glmmTMB objects
##' 
##' Methods have been written that allow \code{glmmTMB} objects to be used with
##' several downstream packages that enable different forms of inference.
##' In particular,
##' \itemize{
##' \item \code{car::Anova} constructs type-II and type-III Anova tables
##' for the fixed effect parameters of the conditional model (this might work with the
##' fixed effects of the zero-inflation or dispersion models, but has not been tested)
##' \item the \code{effects} package computes graphical tabular effect displays
##' (again, for the fixed effects of the conditional component)
##' \item the \code{emmeans} package computes estimated marginal means (aka least-squares means)
##' for the fixed effects of the conditional component
##' }
##' @rdname downstream_methods
##' @param mod a glmmTMB model
##' @param object a glmmTMB model
##' @param trms The \code{terms} component of \code{object} (typically with the
##' response deleted, e.g. via \code{\link{delete.response}}
##' @param xlev Named list of factor levels (\emph{excluding} ones coerced to
##' factors in the model formula)
##' @param \dots Additional parameters that may be supported by the method.
##' @param grid A \code{data.frame} (provided by \code{ref_grid}) containing the
##' predictor settings needed in the reference grid
##' @details While the examples below are disabled for earlier versions of
##' R, they may still work; it may be necessary to refer to private
##' versions of methods, e.g. \code{glmmTMB:::Anova.glmmTMB(model, ...)}.
##' @examples
##' warp.lm <- glmmTMB(breaks ~ wool * tension, data = warpbreaks)
##' if (require(emmeans)) {
##'     emmeans (warp.lm, poly ~ tension | wool)
##' }
##' if (getRversion() >= "3.6.0") {
##'    if (require(car)) {
##'        Anova(warp.lm,type="III")
##'    }
##'    if (require(effects) 
##'        plot(allEffects(warp.lm))
##'    }
##' }


## recover_data method -- DO NOT export -- see zzz.R

#' @importFrom stats delete.response
recover_data.glmmTMB <- function(object, ...) {
    fcall <- getCall(object)
    if (!requireNamespace("emmeans"))
        stop("please install (if necessary) and load the emmeans package")
    emmeans::recover_data(fcall,delete.response(terms(object)),
                 attr(model.frame(object),"na.action"), ...)
}


## emm_basis method -- Dynamically exported, see zzz.R

#' @rdname downstream_methods
#' @aliases downstream_methods
#' @param component which component of the model to compute emmeans for (conditional ("cond"), zero-inflation ("zi"), or dispersion ("disp"))
emm_basis.glmmTMB <- function (object, trms, xlev, grid, component="cond", ...) {
    ## Not needed anymore?
    ## if (component != "cond") warning("only tested for conditional component")
    V <- as.matrix(vcov(object)[[component]])
    misc = list()
    if (family(object)$family=="gaussian") {
        dfargs = list(df = df.residual(object))
        dffun = function(k, dfargs) dfargs$df
    } else {

        dffun = function(k, dfargs) Inf
        dfargs = list()

    }
    fam = switch(component,
                 cond = family(object),
                 zi = list(link = "logit"),
                 disp = list(link = "log"))
    
    misc = emmeans::.std.link.labels(fam, misc)
    ## (used to populate the reminder of response scale)
    contrasts = attr(model.matrix(object), "contrasts")
    ## keep only variables found in conditional fixed effects
    contrasts = contrasts[names(contrasts) %in% all.vars(terms(object))]
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = fixef(object)[[component]]
    if (length(bhat) < ncol(X)) {
        kept = match(names(bhat), dimnames(X)[[2]])
        bhat = NA * X[1, ]
        bhat[kept] = fixef(object)[[component]]
        modmat = model.matrix(trms, model.frame(object), contrasts.arg = contrasts)
        nbasis = estimability::nonest.basis(modmat)
    }  else {
        nbasis = estimability::all.estble
    }
    dfargs = list(df = df.residual(object))
    dffun = function(k, dfargs) dfargs$df
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = dfargs, misc = misc)
}
