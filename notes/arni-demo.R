## source("https://raw.githubusercontent.com/glmmTMB/glmmTMB/master/glmmTMB/R/methods.R")
## source("~/git/glmmTMB/glmmTMB/R/methods.R")

library(lme4)
library(glmmTMB)
data(Owls, package="glmmADMB")
Owls <- transform(Owls, ArrivalTime=c(scale(ArrivalTime,scale=FALSE)),
                  NCalls=SiblingNegotiation)

tmb0 <- glmmTMB(Reaction ~ Days, sleepstudy)
tmb1 <- glmmTMB(Reaction ~ Days + (1|Subject), sleepstudy)
tmb2 <- glmmTMB(Reaction ~ Days + (Days|Subject), sleepstudy)
tmb3 <- glmmTMB(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
sleepstudy$Year <- gl(2, 90, labels=2001:2002)
tmb4 <- glmmTMB(Reaction ~ Days + (1|Subject) + (1|Year), sleepstudy)
tmbz <- glmmTMB(NCalls ~ FoodTreatment + ArrivalTime + (ArrivalTime|Nest) + (1|SexParent),
                ziformula = ~1 + (1|Nest) + (1|SexParent),
                data=Owls, family=poisson(link="log"))

## update(tmb4, .~Days+(1|Subject)+(1|Year))
## update(tmb4, .~Days+(1|Year)+(1|Subject))
## tmb4 <- glmmTMB(Reaction ~ Days + (1|Subject) + (1|Year), sleepstudy)
## tmb4b <- glmmTMB(Reaction ~ Days + (1|Year) + (1|Subject), sleepstudy)

## Z-I with random effects:
## glmmTMB(Reaction ~ Days, ziformula = ~(1|g), sleepstudy)

admb0 <- glmmadmb(Reaction ~ Days, sleepstudy, family="gaussian")
admb1 <- glmmadmb(Reaction ~ Days, sleepstudy, family="gaussian", random=~1|Subject)
admb3 <- glmmadmb(Reaction ~ Days, sleepstudy, family="gaussian",
                  random=~(1|Subject)+(0+Days|Subject))

lme1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy, REML=FALSE)
lme2 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML=FALSE)
lme3 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy, REML=FALSE)
lme4 <- lmer(Reaction ~ Days + (1|Subject) + (1|Year), sleepstudy)

cbpp <- transform(cbpp, prop=incidence/size, obs=factor(1:nrow(cbpp)))
glme.binom <- glmer(prop ~ period + (1|herd) + (1|obs), family=binomial, data=cbpp, weights=size)


## tmb.binom <- glmmTMB(prop ~ period + (1|herd) + (1|obs), family=binomial(), data=cbpp, weights=size)

## fixef.glmmTMB <- function(object, ...)
## {
##   output <- object$obj$env$parList(object$fit$par, object$obj$env$last.par.best)$beta
##   names(output) <- colnames(getME(object,"X"))
##   return(output)
## }

## ranef.glmmTMB <- function(object, se.fit=FALSE, ...)
## {
##   re <- data.frame(object$obj$env$parList()$b)
##   names(re) <- object$modelInfo$reTrms$condList$cnms$Subject
##   rownames(re) <- colnames(getME(object,"Z"))
##   output <- list(re)
##   names(output) <- object$modelInfo$grpVar
##   return(output)
## }

fixef.glmmTMB <- function(object, ...) {
  pl <- object$obj$env$parList(object$fit$par, object$obj$env$last.par.best)
  X <- getME(object, "X")
  ffcond <- structure(pl$beta, names=colnames(X))
  Xzi <- getME(object, "Xzi")
  ffzi <- structure(pl$betazi, names=colnames(Xzi))
  ff <- list(conditional_model=ffcond, zero_inflation=ffzi)
  class(ff) <- "fixef.glmmTMB"
  return(ff)
}

ranef.glmmTMB <- function(object, ...) {
  pl <- object$obj$env$parList(object$fit$par, object$obj$env$last.par.best)
  Z <- getME(object, "Z")
  recond <- structure(pl$b, names=colnames(Z))
  Zzi <- getME(object, "Zzi")
  rezi <- structure(pl$bzi, names=colnames(Zzi))
  re <- list(conditional_model=recond, zero_inflation=rezi)
  class(re) <- "ranef.glmmTMB"
  return(re)
}

vcov.glmmTMB <- function(object, ...)
{
  sdr <- if(!is.null(object$sdr)) object$sdr else update(object,se=TRUE)$sdr
  output <- sdr$cov.fixed
  betas <- which(rownames(output) == "beta")
  output <- output[betas, betas]
  rownames(output) <- colnames(output) <- colnames(getME(object, "X"))
  return(output)
}

## vcov.ranef <- diag(sdr$diag.cov.random)

fixef(tmb0)
fixef(tmb1)
fixef(tmb2)
fixef(tmb3)
fixef(tmb4)
fixef(tmbz)

fixef(lme1)
fixef(lme2)
fixef(lme3)

fixef(admb0)
fixef(admb1)
fixef(admb3)
## Rename to (Intercept) and Days
## betad will be returned by the sigma.glmmTMB() method
## theta will be returned by getME

ranef(tmb0)
ranef(tmb1)
ranef(tmb2)
ranef(tmb3)
ranef(tmb4)
ranef(tmbz)

ranef(lme1)
ranef(lme2)
ranef(lme3)

ranef(admb1)
ranef(admb3)
## ranef() in lme4 returns a list of dataframes, one for each random effect term
## we will format our output in a similar way, but it becomes a list of lists when there is more than 1 non-trivial model
## when ziformula has no random effects in it (e.g. ~0 means no zero inflaction,
##   and ~size means zero inflation is proportional to size on the logit scale) then that is a trivial model for the purpose
##   of random effects

vcov(tmb0)
vcov(tmb1)
vcov(tmb2)
vcov(tmb3)
vcov(tmb4)
vcov(tmbz)

vcov(admb0)
vcov(admb1)
vcov(admb3)

vcov(lme1)
vcov(lme2)
vcov(lme3)
## if the model had zero inflation, then we would be interested in the full vcov
## provide vcov(, full=TRUE) that returns something more than the basic

logLik(tmb0)
logLik(tmb1)
logLik(tmb2)
logLik(tmb3)
logLik(tmb4)
logLik(tmbz)

logLik(lme1)
logLik(lme2)
logLik(lme3)

logLik(admb0)
logLik(admb1)
logLik(admb3)


## the vcov matrix for the random effects should be returned by ranef(model, se.fit=TRUE),
##   i.e., sqrt(model$sdr$diag.cov.random)


tmb1 <- update(tmb1, se=TRUE)
tmb1$sdr$diag.cov.random
sqrt(tmb1$sdr$diag.cov.random)








################################################################################

library(TMB)

input <- glmmTMB(Reaction ~ Days + (1|Subject), sleepstudy, se=TRUE, debug=TRUE)
model <- MakeADFun(input$data, input$parameters, DLL="glmmTMB", random="b",
                   checkParameterOrder=FALSE)
time <- system.time(fit <- nlminb(model$par, model$fn, model$gr))

best <- model$env$last.par.best
rep <- sdreport(model)
