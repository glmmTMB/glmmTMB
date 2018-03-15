## methods for extending emmeans to handle glmmTMB objects

#' @export
recover_data.glmmTMB <- function(object, ...) {
    fcall <- getCall(object)
    recover_data(fcall,delete.response(terms(object)),
                 attr(model.frame(object),"na.action"), ...)
}

#' @export
## FIXME: ignore random effects in 'terms' to avoid warning about
## grouping variables being absent/contrasts ignored
emm_basis.glmmTMB <- function (object, trms, xlev, grid, component="cond", ...) {
    if (component != "cond") warning("only tested for conditional component")
    V <- as.matrix(vcov(object)[[component]])
    dfargs <- NA  ## FIXME: residual df?
    ## use this? misc = .std.link.labels(family(object), misc)
    ## (used to populate the reminder 
    contrasts = attr(model.matrix(object), "contrasts")
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = fixef(object)[[component]]
    if (length(bhat) < ncol(X)) {
        kept = match(names(bhat), dimnames(X)[[2]])
        bhat = NA * X[1, ]
        bhat[kept] = fixef(object)[[component]]
        modmat = model.matrix(trms, model.frame(object), contrasts.arg = contrasts)
        nbasis = estimability::nonest.basis(modmat)
    }  else nbasis = estimability::all.estble
    dfargs = list(df = df.residual(object))
    dffun = function(k, dfargs) dfargs$df

    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = dfargs)
}
