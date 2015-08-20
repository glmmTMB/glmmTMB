
source("https://raw.githubusercontent.com/glmmTMB/glmmTMB/master/glmmTMB/R/methods.R")

fixef.glmmTMB <- function(object, ...)
{
  output <- object$obj$env$parList()$beta
  names(output) <- colnames(getME(object,"X"))
  return(output)
}

ranef.glmmTMB <- function(object, se.fit, ...)
{
  output <- object$obj$env$parList()$b
  names(output) <- # we need 'grpVar' from glmmTMB:::getReStruc()
  return(output)
}

vcov.glmmTMB <- function(object, ...)
{
  sdr <- if(!is.null(object$sdr)) object$sdr else update(object,se=TRUE)$sdr
  output <- sdr$cov.fixed
  return(output)
}

## vcov.ranef <- diag(sdr$diag.cov.random)

output <- glmmTMB(Reaction ~ Days + (1|Subject), sleepstudy, se=TRUE)

fixef(output)
## Rename to (Intercept) and Days
## betad will be returned by the sigma.glmmTMB() method
## theta will be returned by getME

ranef(output)
## ranef() in lme4 returns a list of dataframes, one for each random effect term
## we will format our output in a similar way, but it becomes a list of lists when there is more than 1 non-trivial model
## when ziformula has no random effects in it (e.g. ~0 means no zero inflaction,
##   and ~size means zero inflation is proportional to size on the logit scale) then that is a trivial model for the purpose
##   of random effects

vcov(output)[[1]][1:2,1:2]
## if the model had zero inflation, then we would be interested in the full vcov
## provide vcov(, full=TRUE) that returns something more than the basic

round(vcov(output)[[2]])
## the vcov matrix for the random effects should be returned by ranef(model, se.fit=TRUE),
##   i.e., sqrt(model$sdr$diag.cov.random)


model.2 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
class(ranef(model.2)$Subject)
output$sdr$diag.cov.random
sqrt(output$sdr$diag.cov.random)








################################################################################

library(TMB)

input <- glmmTMB(Reaction ~ Days + (1|Subject), sleepstudy, se=TRUE, debug=TRUE)
model <- MakeADFun(input$data, input$parameters, DLL="glmmTMB", random="b",
                   checkParameterOrder=FALSE)
time <- system.time(fit <- nlminb(model$par, model$fn, model$gr))

best <- model$env$last.par.best
rep <- sdreport(model)
