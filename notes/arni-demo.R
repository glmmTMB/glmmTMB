## source("https://raw.githubusercontent.com/glmmTMB/glmmTMB/master/glmmTMB/R/methods.R")
source("~/git/glmmTMB/glmmTMB/R/methods.R")

## library(glmmTMB)
library(reshape)
data(Owls, package="glmmADMB")
Owls <- rename(Owls, c(SiblingNegotiation="NCalls"))
Owls <- transform(Owls, ArrivalTime=scale(ArrivalTime,center=TRUE,scale=FALSE))

tmb0 <- glmmTMB(Reaction ~ Days, sleepstudy)
tmb1 <- glmmTMB(Reaction ~ Days + (1|Subject), sleepstudy)
tmb2 <- glmmTMB(Reaction ~ Days + (Days|Subject), sleepstudy)
tmb3 <- glmmTMB(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
tmbz <- glmmTMB(NCalls ~ FoodTreatment * SexParent + ArrivalTime * SexParent
                + offset(logBroodSize) + (1|Nest), ziformula = ~1,
                data=Owls, family=poisson(link="log"))

## Z-I with random effects:
## glmmTMB(Reaction ~ Days, ziformula = ~(1|g), sleepstudy)

admb0 <- glmmadmb(Reaction ~ Days, sleepstudy, family="gaussian")
admb1 <- glmmadmb(Reaction ~ Days, sleepstudy, family="gaussian", random=~1|Subject)
admb3 <- glmmadmb(Reaction ~ Days, sleepstudy, family="gaussian", random=~(1|Subject)+(0+Days|Subject))

lme1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy, REML=FALSE)
lme2 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML=FALSE)
lme3 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy, REML=FALSE)

## fixef.glmmTMB <- function(object, ...)
## {
##   output <- object$obj$env$parList(object$fit$par, object$obj$env$last.par.best)$beta
##   names(output) <- colnames(getME(object,"X"))
##   return(output)
## }

ranef.glmmTMB <- function(object, se.fit=FALSE, ...)
{
  re <- data.frame(object$obj$env$parList()$b)
  names(re) <- object$modelInfo$reTrms$condList$cnms$Subject
  rownames(re) <- colnames(getME(object,"Z"))
  output <- list(re)
  names(output) <- object$modelInfo$grpVar
  return(output)
}

vcov.glmmTMB <- function(object, ...)
{
  sdr <- if(!is.null(object$sdr)) object$sdr else update(object,se=TRUE)$sdr
  output <- sdr$cov.fixed
  betas <- which(rownames(output) == "beta")
  output <- output[betas, betas]
  rownames(output) <- colnames(output) <- colnames(getME(object,"X"))
  return(output)
}

## vcov.ranef <- diag(sdr$diag.cov.random)

fixef(tmb0)
fixef(tmb1)
fixef(tmb2)
fixef(tmb3)

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

ranef(lme1)
ranef(lme2)
ranef(lme3)

ranef(admb0)
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

vcov(admb0)
vcov(admb1)
vcov(admb3)

vcov(lme1)
vcov(lme2)
vcov(lme3)
## if the model had zero inflation, then we would be interested in the full vcov
## provide vcov(, full=TRUE) that returns something more than the basic

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
