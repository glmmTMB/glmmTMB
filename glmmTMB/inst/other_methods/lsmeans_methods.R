## EXPERIMENTAL (not yet exported)
#' @importFrom stats delete.response model.frame na.pass
recover.data.glmmTMB <- function(object, ...) {
    fcall <- getCall(object)
    lsmeans::recover.data(fcall,delete.response(terms(object)),
                 attr(model.frame(object),"na.action"), ...)
}
lsm.basis.glmmTMB <- function (object, trms, xlev, grid, vcov.,
                               mode = "asymptotic", component="cond", ...) {
    if (mode != "asymptotic") stop("only asymptotic mode is available")
    if (component != "cond") stop("only tested for conditional component")
    if (missing(vcov.)) 
        V <- as.matrix(vcov(object)[[component]])
    else V <- as.matrix(.my.vcov(object, vcov.))
    dfargs = misc = list()
    if (mode == "asymptotic") {
        dffun = function(k, dfargs) NA
    }
    ## use this? misc = .std.link.labels(family(object), misc)
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
    }
    else nbasis = estimability::all.estble
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
        dfargs = dfargs, misc = misc)
}
