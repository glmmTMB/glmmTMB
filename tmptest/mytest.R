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

predict.glmmTMB <- function(object,newdata=NULL,...) {
   
  if (is.null(newdata)) {
    stop("newdata=NULL case not yet written")
    ## FIXME: in sdr object
  }
  mf <- object$call
  ## FIXME: DRY so much
  ## now work on evaluating model frame
  ## do we want to re-do this part???
  m <- match(c("subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$data <- newdata
  mf[[1]] <- as.name("model.frame")

  mf$formula <- object$modelInfo$allForm$combForm
  newFr <- eval(mf, parent.frame())
  respCol <- match(respNm <- names(object$modelInfo$respCol),names(newFr))
  ## create *or* overwrite response column for prediction data with NA
  newFr$respNm <- NA
  
  ## FIXME: not yet handling population-level predictions (re.form
  ##  or new levels/allow.new.levels)
  
  ## append to existing model frame
  augFr <- rbind(object$fr,newFr)
  
  ## now re-do 
}