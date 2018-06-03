## modified from car::effectsmer.R

## mer.to.glm evaluates a 'glm' model that is as similar to a given 'mer'
## model as possible.  It is of class c("fakeglm", "glm", "lm")
## several items are added to the created objects. Do not export

#' @importFrom lme4 nobars
## glmmTMB.to.glm <- function(mod, KR=FALSE) {
##     if (KR) { ##  && !requireNamespace("pbkrtest", quietly=TRUE)){
##         KR <- FALSE
##         warning("pbkrtest is not compatible with glmmTMB, KR set to FALSE")
##     }
##     orig_family <- family(mod)
##     link <- orig_family$link
##     family <- orig_family$family
##     cl <- getCall(mod)
##     ## prevents evaluation (suggested by Nate TeGrotenhuis)
##     cl$control <- glm.control(epsilon=1) 
##     .m <- match(c("formula", "family", "data", "weights", "subset",
##                   "na.action", "offset",
##                   "model", "contrasts"), names(cl), 0L)
##     cl <- cl[c(1L, .m)]
##     cl[[1L]] <- as.name("glm")
##     cl$formula <- lme4::nobars(as.formula(cl$formula))
##     ##    cl$data <- mod@frame # caused bug with a 'poly' in the formula
##     cl$family <- gaussian
##     mod2 <- eval(cl)
##     cl$family <- orig_family
##     mod2$coefficients <- glmmTMB::fixef(mod)[["cond"]]
##     ## mod2$vcov <- if (family == "gaussian" && link == "identity" && KR) as.matrix(pbkrtest::vcovAdj(mod)) else as.matrix(vcov(mod))
##     mod2$vcov <- vcov(mod)[["cond"]]
##     mod2$linear.predictors <- model.matrix(mod2) %*% mod2$coefficients
##     mod2$fitted.values <- mod2$family$linkinv(mod2$linear.predictors)
##     mod2$weights <- as.vector(with(mod2,
##                                    prior.weights * (family$mu.eta(linear.predictors)^2 /
##                                                     family$variance(fitted.values))))
##     mod2$residuals <- with(mod2,
##                            prior.weights * (y - fitted.values)/weights )
##     class(mod2) <- c("fakeglm", class(mod2))
##     return(mod2)
## }

##method for 'fakeglm' objects. Do not export
vcov.fakeglm <- function(object, ...) object$vcov

##The next six functions should be exported as S3 methods

#' @export
Effect.glmmTMB <- function (focal.predictors, mod, ...,
                            component="cond", KR = FALSE)  {
    if (KR && !requireNamespace("pbkrtest", quietly = TRUE)) {
        KR <- FALSE
        warning("pbkrtest is not available, KR set to FALSE")
    }
    fam <- family(mod)
    args <- list(call = getCall(mod),
                 coefficients = lme4::fixef(mod)[[component]],
                 vcov = vcov(mod)[[component]],
                 family=family(mod))
    Effect.default(focal.predictors, mod, ..., sources = args)
}

## testing
if (FALSE) {
    ## working directory: inst/other_methods/
    ## source("effectsglmmTMB.R")
    library(effects)
    gg1 <- glmmTMB(round(Reaction)~Days+(1|Subject),family=poisson,data=lme4::sleepstudy)
    ## glmmTMB.to.glm(gg1)
    allEffects(gg1)
    set.seed(101)
    dd <- data.frame(y=rnbinom(1000,mu=4,size=1),
                     x = rnorm(1000),
                     f=factor(rep(LETTERS[1:20],each=50)))
    gg2 <- glmmTMB(y~x+(1|f),family=nbinom2,data=dd)
    ls(environment(gg2$modelInfo$family$dev.resids),all.names=TRUE)
    ## debug(gg2$modelInfo$family$dev.resids)
    allEffects(gg2)
}
