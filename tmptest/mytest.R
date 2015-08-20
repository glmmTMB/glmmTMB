## Demonstrate principles behind predictions

library(glmmTMB)
data(sleepstudy,package="lme4")    
qw <- glmmTMB(Reaction~Days+(Days|Subject),sleepstudy)

## Simple predict method:
library(TMB) ## FIXME!
H <- optimHess(qw$fit$par,qw$obj$fn,qw$obj$gr)
sdr <- sdreport(qw$obj, hessian=H)  ## FIXME: Possible to ignore.parm.uncertainty=TRUE
summary(sdr,"report")

## Add newdata
options(na.action=na.pass)
newdata <- data.frame(Reaction = NA, Days = 0:9, Subject = "373")
sleepstudy2 <- rbind( sleepstudy, newdata )
qw2 <- glmmTMB(Reaction~Days+(Days|Subject),sleepstudy2) ## FIXME: NO REFIT !!!
sdr2 <- sdreport(qw2$obj, par=qw$fit$par, hessian=H) ## NOTE: par + hessian from previous fit
summary(sdr2,"report")

head(summary(sdr2,"report"),10)
tail(summary(sdr2,"report"),10)
