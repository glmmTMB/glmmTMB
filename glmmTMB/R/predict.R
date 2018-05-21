## Helper function for predict.
## Assert that we can use old model (data.tmb0) as basis for
## predictions using the new data (data.tmb1):
assertIdenticalModels <- function(data.tmb1, data.tmb0, allow.new.levels=FALSE)
{
    ## Check terms. Only 'blockReps' and 'blockSize' are allowed to
    ## change.  Note that we allow e.g. spatial covariance matrices to
    ## change, while e.g. an unstrucured covariance must remain the
    ## same.
    checkTerms <- function(t1, t0) {
        ## Defensive check:
        stopifnot(identical(names(t1), names(t0)))
        ## *Never* allowed to differ:
        testIdentical <- function(checkNm) {
            unlist( Map( function(x,y)
                identical(x[checkNm], y[checkNm]), t0, t1) )
        }
        ok <- testIdentical( c("blockNumTheta", "blockCode") )
        if ( ! all(ok) ) {
            msg <- c("Prediction is not possible for terms: ",
                     paste(names(t1)[!ok], collapse=", "), "\n",
                     "Probably some factor levels in 'newdata' require fitting a new model.")
            stop(msg)
        }
        ## Sometimes allowed to differ:
        if ( ! allow.new.levels ) {
            ok <- testIdentical( c( "blockReps", "blockSize") )
            if ( ! all(ok) ) {
                msg <- c("Predicting new random effect levels for terms: ",
                         paste(names(t1)[!ok], collapse=", "), "\n",
                         "Disable this warning with 'allow.new.levels=TRUE'")
                ## FIXME: warning or error ?
                warning(msg)
            }
        }
    }
    checkTerms( data.tmb1$terms,   data.tmb0$terms )
    checkTerms( data.tmb1$termszi, data.tmb0$termszi )
    ## Fixed effect parameters must be identical
    checkModelMatrix <- function(X1, X0) {
        if( !identical(colnames(X1), colnames(X0)) ) {
            msg <- c("Prediction is not possible for unknown fixed effects: ",
                     paste( setdiff(colnames(X1), colnames(X0)), collapse=", "), "\n",
                     "Probably some factor levels in 'newdata' require fitting a new model.")
            stop(msg)
        }
    }
    checkModelMatrix(data.tmb1$X,   data.tmb0$X)
    checkModelMatrix(data.tmb1$Xzi, data.tmb0$Xzi)
    NULL
}

##' prediction
##' @param object a \code{glmmTMB} object
##' @param newdata new data for prediction
##' @param se.fit return the standard errors of the predicted values?
##' @param zitype deprecated: formerly used to specify type of zero-inflation probability. Now synonymous with \code{type}
##' @param type Denoting \eqn{mu} as the mean of the conditional distribution and
##' \code{p} as the zero-inflation probability,
##' the possible choices are:
##' \describe{
##' \item{"link"}{conditional mean on the scale of the link function,
##' or equivalently the linear predictor of the conditional model}
##' \item{"response"}{expected value; this is \eqn{mu*(1-p)} for zero-inflated models
##' and \code{mu} otherwise}
##' \item{"conditional"}{mean of the conditional response; \code{mu} for all models
##' (i.e., synonymous with \code{"response"} in the absence of zero-inflation}
##' \item{"zprob"}{the probability of a structural zero (gives an error
##' for non-zero-inflated models)}
##' ##' \item{"zilink"}{predicted zero-inflation probability on the scale of
##' the logit link function}
##' }
##' @param na.action how to handle missing values in \code{newdata} (see \code{\link{na.action}});
##' the default (\code{na.pass}) is to predict \code{NA}
##' @param debug (logical) return the \code{TMBStruc} object that will be
##' used internally for debugging?
##' @param re.form (not yet implemented) specify which random effects to condition on when predicting
##' @param allow.new.levels allow previously unobserved levels in random-effects variables? see details.
##' @param \dots unused - for method compatibility
##' @details
##' Prediction of new random effect levels is possible as long as the model specification (fixed effects and parameters) is kept constant.
##' However, to ensure intentional usage, a warning is triggered if \code{allow.new.levels=FALSE} (the default).
##' @examples
##' data(sleepstudy,package="lme4")
##' g0 <- glmmTMB(Reaction~Days+(Days|Subject),sleepstudy)
##' predict(g0, sleepstudy)
##' ## Predict new Subject
##' nd <- sleepstudy[1,]
##' nd$Subject <- "new"
##' predict(g0, newdata=nd, allow.new.levels=TRUE)
##' @importFrom TMB sdreport
##' @importFrom stats optimHess model.frame na.fail na.pass napredict
##' @export
predict.glmmTMB <- function(object,newdata=NULL,
                            se.fit=FALSE,
                            re.form, allow.new.levels=FALSE,
                            type = c("link", "response",
                                     "conditional","zprob","zlink"),
                            zitype = NULL,
                            na.action = na.pass,
                            debug=FALSE,
                            ...)
{
  ## FIXME: add re.form

  if (!is.null(zitype)) {
     warning("zitype is deprecated: please use type instead")
     type <- zitype
  }
  type <- match.arg(type)
  if (!missing(re.form)) stop("re.form not yet implemented")
  ##if (allow.new.levels) stop("allow.new.levels not yet implemented")
  mc <- mf <- object$call
  ## FIXME: DRY so much
  ## now work on evaluating model frame
  ## do we want to re-do this part???

  ## need to 'fix' call to proper model.frame call whether or not
  ## we have new data, because ... (??)
  m <- match(c("subset", "weights", "offset", "na.action"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]

  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf$formula <- RHSForm(object$modelInfo$allForm$combForm, as.form=TRUE)
    
  if (is.null(newdata)) {
    mf$data <- mc$data ## restore original data
    newFr <- object$fr
  } else {
    mf$data <- newdata
    mf$na.action <- na.action
    newFr <- eval.parent(mf)
  }

  omi <- object$modelInfo  ## shorthand

  respCol <- match(respNm <- names(omi$respCol),names(newFr))
  ## create *or* overwrite response column for prediction data with NA
  newFr[[respNm]] <- NA

  ## FIXME: not yet handling population-level predictions (re.form
  ##  or new levels/allow.new.levels)

  ## append to existing model frame
  augFr <- rbind(object$fr,newFr)

  ## Pointers into 'new rows' of augmented data frame.
  w <- nrow(object$fr) + seq_len(nrow(newFr))

  ## Variety of possible binomial inputs are taken care of by
  ## 'mkTMBStruc' further down.
  yobs <- augFr[[names(omi$respCol)]]

  ## match type arg with internal name
  ## FIXME: warn if "link"  
  ziPredNm <- switch(type,
                     response   = "corrected",
                     link       =,
                     conditional= "uncorrected",
                     zlink      = ,
                     zprob      = "prob",
                     stop("unknown type ",type))
  ziPredCode <- .valid_zipredictcode[ziPredNm]

  ## need eval.parent() because we will do eval(mf) down below ...
  TMBStruc <-
        ## FIXME: make first arg of mkTMBStruc into a formula list
        ## with() interfering with eval.parent() ?
        eval.parent(mkTMBStruc(RHSForm(omi$allForm$formula,as.form=TRUE),
                               omi$allForm$ziformula,
                               omi$allForm$dispformula,
                               omi$allForm$combForm,
                               mf,
                               fr=augFr,
                               yobs=yobs,
                               respCol=respCol,
                               weights=model.weights(augFr),
                               family=omi$family,
                               ziPredictCode=ziPredNm,
                               doPredict=as.integer(se.fit),
                               whichPredict=w))

  ## short-circuit
  if(debug) return(TMBStruc)

  ## Check that the model specification is unchanged:
  assertIdenticalModels(TMBStruc$data.tmb,
                        object$obj$env$data, allow.new.levels)

  newObj <- with(TMBStruc,
                 MakeADFun(data.tmb,
                           parameters,
                           map = mapArg,
                           random = randomArg,
                           profile = NULL, # TODO: Optionally "beta"
                           silent = TRUE,
                           DLL = "glmmTMB"))

  oldPar <- object$fit$par
  newObj$fn(oldPar)  ## call once to update internal structures
  lp <- newObj$env$last.par

  na.act <- attr(model.frame(object),"na.action")
  do.napred <- missing(newdata) && !is.null(na.act)
  if (!se.fit) {
      pred <- newObj$report(lp)$mu_predict
  } else {
      H <- with(object,optimHess(oldPar,obj$fn,obj$gr))
      ## FIXME: Eventually add 'getReportCovariance=FALSE' to this sdreport
      ##        call to fix memory issue (requires recent TMB version)
      ## Fixed! (but do we want a flag to get it ? ...)
      sdr <- sdreport(newObj,oldPar,hessian.fixed=H,getReportCovariance=FALSE)
      sdrsum <- summary(sdr, "report") ## TMB:::summary.sdreport(sdr, "report")
      pred <- sdrsum[,"Estimate"]
      se <- sdrsum[,"Std. Error"]
  }
  if (do.napred) {
      pred <- napredict(na.act,pred)
      if (se.fit) se <- napredict(na.act,se)
  }
  if (type %in% c("zlink","link")) {
     ff <- object$modelInfo$family
     if (!(type=="link" && ff$link=="identity")) {
         if (type=="zlink") {
             ff <- make.link("logit")
         }
         pred <- ff$linkfun(pred)
         if (se.fit) se <- se/ff$mu.eta(pred) ## do this after transforming pred!
     } ## if not identity link  
  } ## if link or zlink
  if (!se.fit) return(pred) else return(list(fit=pred,se.fit=se))
}
