##' prediction
##' @param object a \code{glmmTMB} object
##' @param newdata new data for prediction
##' @param se.fit return the standard errors of the predicted values?
##' @param zitype for zero-inflated models,
##' return expected value ("response": (mu*(1-p))),
##' the mean of the conditional distribution ("conditional": mu),
##' or the probability of a structural zero ("zprob")?
##' @param debug (logical) return the \code{TMBStruc} object that will be
##' used internally for debuggin?
##' @param re.form (not yet implemented) specify which random effects to condition on when predicting
##' @param allow.new.levels (not yet implemented) allow previously unobserved levels in random-effects grouping variables?
##' @param \dots unused - for method compatibility

##' @examples
##' data(sleepstudy,package="lme4")
##' g0 <- glmmTMB(Reaction~Days+(Days|Subject),sleepstudy)
##' predict(g0, sleepstudy)
##' @importFrom TMB sdreport
##' @importFrom stats optimHess
##' @export
predict.glmmTMB <- function(object,newdata=NULL,
                            se.fit=FALSE,
                            re.form, allow.new.levels=FALSE,
                            zitype = c("response","conditional","zprob"),
                            debug=FALSE,
                            ...)
{
  ## FIXME: add re.form, type, ...
  ## FIXME: deal with napredict stuff ...

  if (!missing(re.form)) stop("re.form not yet implemented")
  if (allow.new.levels) stop("allow.new.levels not yet implemented")
  mf <- object$call
  ## FIXME: DRY so much
  ## now work on evaluating model frame
  ## do we want to re-do this part???

  if (is.null(newdata)) {
      newFr <- object$fr
  } else {
      m <- match(c("subset", "weights", "na.action", "offset"),
                 names(mf), 0L)
      mf <- mf[c(1L, m)]
      mf$drop.unused.levels <- TRUE
      mf$data <- newdata
      mf[[1]] <- as.name("model.frame")

      mf$formula <- object$modelInfo$allForm$combForm
      newFr <- eval(mf, parent.frame())
  }
  respCol <- match(respNm <- names(object$modelInfo$respCol),names(newFr))
  ## create *or* overwrite response column for prediction data with NA
  newFr[[respNm]] <- NA

  ## FIXME: not yet handling population-level predictions (re.form
  ##  or new levels/allow.new.levels)

  ## append to existing model frame
  augFr <- rbind(object$fr,newFr)

  yobs <- augFr[[names(object$modelInfo$respCol)]]

  ## match zitype arg with internal name
  ziPredNm <- switch(match.arg(zitype),
                       response="corrected",
                       conditional="uncorrected",
                         zprob="prob",
                       stop("unknown zitype ",zitype))
  ziPredCode <- .valid_zipredictcode[ziPredNm]
  TMBStruc <- with(object$modelInfo,
                   ## FIXME: make first arg of mkTMBStruc into a formula list
                   mkTMBStruc(allForm$formula,
                              allForm$ziformula,allForm$dispformula,
                         mf,augFr,
                         yobs=augFr[[names(respCol)]],
                         offset=NULL,weights=NULL,
                         family=familyStr,link=link,
                         ziPredictCode=ziPredNm,
                         doPredict=1))

  ## short-circuit
  if(debug) return(TMBStruc)

  newObj <- with(TMBStruc,
                 MakeADFun(data.tmb,
                           parameters,
                           random = randomArg,
                           profile = NULL, # TODO: Optionally "beta"
                           silent = TRUE,
                           DLL = "glmmTMB"))

  oldPar <- object$fit$par
  newObj$fn(oldPar)  ## call once to update internal structures
  H <- with(object,optimHess(oldPar,obj$fn,obj$gr))
  sdr <- sdreport(newObj,oldPar,hessian.fixed=H)
  sdr.rpt <- summary(sdr, "report") ## TMB:::summary.sdreport(sdr, "report")
  ## now strip off original values
  w <- which(is.na(augFr[[respNm]]))

  pred <- sdr.rpt[w, , drop=FALSE]

  if (!se.fit) {
      return(pred[,"Estimate"])
  } else {
      return(list(fit=pred[,"Estimate"],
                  se.fit=pred[,"Std. Error"]))
  }
}
