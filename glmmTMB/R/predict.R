##' prediction
##' @param object a \code{glmmTMB} object
##' @param newdata new data for prediction
##' @examples 
##' data(sleepstudy,package="lme4")
##' g0 <- glmmTMB(Reaction~Days+(Days|Subject),sleepstudy)
##' @export
predict.glmmTMB <- function(object,newdata=NULL,...) {
  ## FIXME: add re.form, type, ...
  
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