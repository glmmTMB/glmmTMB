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
##' 
##' In particular,
##' \itemize{
##' \item \code{car::Anova} constructs type-II and type-III Anova tables
##' for the fixed effect parameters of any component
##' \item the \code{emmeans} package computes estimated marginal means (previously known as least-squares means)
##' for the fixed effects of any component
##' \item the \code{effects} package computes graphical tabular effect displays
##' (only for the fixed effects of the conditional component)
##' }
##' @param mod a glmmTMB model
##' @param component which component of the model to test/analyze ("cond", "zi", or "disp")
##' @param \dots Additional parameters that may be supported by the method.
##' @details While the examples below are disabled for earlier versions of
##' R, they may still work; it may be necessary to refer to private
##' versions of methods, e.g. \code{glmmTMB:::Anova.glmmTMB(model, ...)}.
##' @importFrom stats delete.response
##' @examples
##' warp.lm <- glmmTMB(breaks ~ wool * tension, data = warpbreaks)
##' salamander1 <- up2date(readRDS(system.file("example_files","salamander1.rds",package="glmmTMB")))
##' if (require(emmeans)) {
##'     emmeans(warp.lm, poly ~ tension | wool)
##'     emmeans(salamander1, ~ mined, type="response")
##'     emmeans(salamander1, ~ mined, component="zi", type="response")
##' }
##' if (getRversion() >= "3.6.0") {
##'    if (require(car)) {
##'        Anova(warp.lm,type="III")
##'        Anova(salamander1)
##'        Anova(salamander1, component="zi")
##'    }
##'    if (require(effects)) {
##'        plot(allEffects(warp.lm))
##'        plot(allEffects(salamander1))
##'    }
##' }
NULL  ## don't document the files here!


## recover_data method -- DO NOT export -- see zzz.R
## do not document either

recover_data.glmmTMB <- function(object, ...) {
    fcall <- getCall(object)
    if (!requireNamespace("emmeans"))
        stop("please install (if necessary) and load the emmeans package")
    emmeans::recover_data(fcall,delete.response(terms(object)),
                 attr(model.frame(object),"na.action"), ...)
}



## emm_basis method -- Dynamically exported, see zzz.R
## don't document, causes confusion

## @rdname downstream_methods
## @aliases downstream_methods
## @param component which component of the model to compute emmeans for (conditional ("cond"), zero-inflation ("zi"), or dispersion ("disp"))
## vcov. user-specified covariance matrix
emm_basis.glmmTMB <- function (object, trms, xlev, grid, component="cond", vcov., ...) {
    ## browser()
    L <- list(...)
    if (length(L)>0) {
        ## don't warn: $misc and $options are always passed through ...
        ## warning("ignored extra arguments to emm_basis.glmmTMB: ",
        ## paste(names(L),collapse=", "))
    }
    if (missing(vcov.)) {
        V <- as.matrix(vcov(object)[[component]])
    } else {
        V <- vcov.
    }
    misc <- list()
    if (family(object)$family=="gaussian") {
        dfargs = list(df = df.residual(object))
        dffun = function(k, dfargs) dfargs$df
    } else {

        dffun = function(k, dfargs) Inf
        dfargs = list()

    }
    fam <- switch(component,
                 cond = family(object),
                 zi = list(link="logit"),
                 disp = list(link="log"))
    
    misc <- emmeans::.std.link.labels(fam, misc)
    ## (used to populate the reminder of response scale)
    contrasts <- attr(model.matrix(object), "contrasts")
    ## keep only variables found in conditional fixed effects
    contrasts <- contrasts[names(contrasts) %in% all.vars(terms(object))]
    m <- model.frame(trms, grid, na.action=na.pass, xlev=xlev)
    X <- model.matrix(trms, m, contrasts.arg=contrasts)
    bhat <- fixef(object)[[component]]
    if (length(bhat) < ncol(X)) {
        kept <- match(names(bhat), dimnames(X)[[2]])
        bhat <- NA * X[1, ]
        bhat[kept] <- fixef(object)[[component]]
        modmat <- model.matrix(trms, model.frame(object), contrasts.arg=contrasts)
        nbasis <- estimability::nonest.basis(modmat)
    }  else {
        nbasis <- estimability::all.estble
    }
    dfargs <- list(df=df.residual(object))
    dffun <- function(k, dfargs) dfargs$df
    namedList(X, bhat, nbasis, V, dffun, dfargs, misc)
}
