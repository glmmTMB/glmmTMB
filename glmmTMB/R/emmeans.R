## methods for extending emmeans to handle glmmTMB objects

#' @rdname emmeans.glmmTMB
#' @examples
#' if (require(emmeans)) {
#'    warp.lm <- glmmTMB(breaks ~ wool * tension, data = warpbreaks)
#'    emmeans (warp.lm, poly ~ tension | wool)
#' }
#' @export
recover_data.glmmTMB <- function(object, ...) {
    fcall <- getCall(object)
    recover_data(fcall,delete.response(terms(object)),
                 attr(model.frame(object),"na.action"), ...)
}

## copied from emmeans (not exported)
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

#' @rdname emmeans.glmmTMB
#' @export
## FIXME: ignore random effects in 'terms' to avoid warning about
## grouping variables being absent/contrasts ignored
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
    
    ## use this? misc = .std.link.labels(family(object), misc)
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
    }  else nbasis = estimability::all.estble
    dfargs = list(df = df.residual(object))
    dffun = function(k, dfargs) dfargs$df

    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = dfargs)
}
