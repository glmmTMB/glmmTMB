## modified from contribution by Sanford Weisberg

##' @rdname downstream_methods
##' @param focal.predictors a character vector of one or more predictors in the
##'  model in any order.

##'
##' @rawNamespace if(getRversion() >= "3.6.0") {
##'   S3method(effects::Effect, glmmTMB)
##' } else {
##'   export(Effect.glmmTMB)
##' }
Effect.glmmTMB <- function (focal.predictors, mod, ...) {
    fam <- family(mod)
    ## code to make the 'truncated_*' families work
    if (grepl("^truncated", fam$family)) 
        fam <- c(fam, make.link(fam$link))
    ## dummy functions to make Effect.default work
    dummyfuns <- list(variance=function(mu) mu,
                      initialize=expression(mustart <- y + 0.1),
                      dev.resids=function(...) poisson()$dev.res(...)
                      )
    for (i in names(dummyfuns)) {
        if (is.null(fam[[i]])) fam[[i]] <- dummyfuns[[i]]
    }
    ## allow calculation of effects ...
    if (length(formals(fam$variance))>1) {
        warning("overriding variance function for effects: ",
                "computed variances may be incorrect")
        fam$variance <- dummyfuns$variance
    }
    args <- list(call = getCall(mod),
                 coefficients = lme4::fixef(mod)[["cond"]],
                 vcov = vcov(mod)[["cond"]],
                 family=fam)
    if (!requireNamespace("effects"))
        stop("please install the effects package")
    effects::Effect.default(focal.predictors, mod, ..., sources = args)
}
