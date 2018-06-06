## modified from car::effectsmer.R

globalVariables("Effect") ## suppress code check warning

##' @rdname downstream_methods
##' @param focal.predictors a character vector of one or more predictors in the
##'  model in any order.

##' @export Effect.glmmTMB
Effect.glmmTMB <- function (focal.predictors, mod, ...,
                            component="cond")  {
    fam <- family(mod)
    args <- list(call = getCall(mod),
                 coefficients = lme4::fixef(mod)[[component]],
                 vcov = vcov(mod)[[component]],
                 family=family(mod))
    if (!requireNamespace("effects"))
        stop("please install the effects package")
    ## use unclass() to get to Effect.default
    Effect(focal.predictors, unclass(mod), ..., sources = args)

}

