## methods for extending emmeans to handle glmmTMB objects

## FIXME:
## I don't know a good way to export methods
## conditionally defining a generic as described in this thread
## <https://stat.ethz.ch/pipermail/r-package-devel/2018q2/002764.html>
## seems like a good idea, but loading the other package masks the
## generic!
## so: fully export methods

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
##' @examples
##' warp.lm <- glmmTMB(breaks ~ wool * tension, data = warpbreaks)
##' if (require(emmeans)) {
##'     emmeans (warp.lm, poly ~ tension | wool)
##' }
##' if (require(car)) {
##'     Anova(warp.lm,type="III")
##' }
##' if (require(effects)) {
##'     plot(allEffects(warp.lm))
##' }


## recover_data <- function(mod, ...) {
##    if (requireNamespace("emmeans", quietly = TRUE)) {
##        emmeans::recover_data(object, ...)
##    } else UseMethod("recover_data")
## }

## ##' rdname downstream_methods
## ##' export 
## emm_basis <-  function(object, trms, xlev, grid, ...) {
##     if (requireNamespace("emmeans", quietly = TRUE)) {
##         emmeans::emm_basis(object, trms, xlev, grid, ...)
##     } else UseMethod("emm_basis")
## }

#' @importFrom stats delete.response
#' @rawNamespace if(getRversion() >= "3.6.0") {
#'   S3method(emmeans::recover_data, glmmTMB)
#' } else {
#'   export(recover_data.glmmTMB)
#' }
recover_data.glmmTMB <- function(object, ...) {
    fcall <- getCall(object)
    if (!requireNamespace("emmeans"))
        stop("please install (if necessary) and load the emmeans package")
    emmeans::recover_data(fcall,delete.response(terms(object)),
                 attr(model.frame(object),"na.action"), ...)
}

## copied from emmeans (not exported)
## (will be exported in next release of emmeans)
.std.link.labels <- function (fam, misc) 
{
    if (is.null(fam) || !is.list(fam)) 
        return(misc)
    if (fam$link == "identity") 
        return(misc)
    misc$tran = fam$link
    misc$inv.lbl = "response"
    if (grepl("binomial", fam$family))  {
        misc$inv.lbl = "prob"
    } else if (fam$family=="beta") {
        misc$inv.lbl = "prop"
    } else if (grepl("(pois|nbinom|tweedie|Gamma)", fam$family)) {
        misc$inv.lbl = "rate"
    }
    misc
}

#' @rdname downstream_methods
#' @aliases downstream_methods
#' @param component which component of the model to compute emmeans for (conditional ("cond"), zero-inflation ("zi"), or dispersion ("disp"))

#' @rawNamespace if(getRversion() >= "3.6.0") {
#'   S3method(emmeans::emm_basis, glmmTMB)
#' } else {
#'   export(emm_basis.glmmTMB)
#' }
emm_basis.glmmTMB <- function (object, trms, xlev, grid, component="cond", ...) {
    if (component != "cond") warning("only tested for conditional component")
    V <- as.matrix(vcov(object)[[component]])
    misc = list()
    if (family(object)$family=="gaussian") {
        dfargs = list(df = df.residual(object))
        dffun = function(k, dfargs) dfargs$df
    } else {

        dffun = function(k, dfargs) Inf
        dfargs = list()

    }
    
    misc = .std.link.labels(family(object), misc)
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
