##' Extract the fixed-effects estimates
##'
##' Extract the estimates of the fixed-effects parameters from a fitted model.
##' @name fixef
##' @title Extract fixed-effects estimates
##' @aliases fixef fixed.effects fixef.merMod
##' @docType methods
##' @param object any fitted model object from which fixed effects estimates can
##' be extracted.
##' @param \dots optional additional arguments. Currently none are used in any
##' methods.
##' @return a named, numeric vector of fixed-effects estimates.
##' @keywords models
##' @examples
##' fixef(glmmtmb(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
##' @importFrom nlme fixef
##' @export fixef
##' @method fixef glmmTMB
##' @export
fixef.glmmTMB <- function(object,...) {
  X <- getME(object,"X")
  pl <- object$obj$env$parList(object$fit$par, object$obj$env$last.par.best)
  ffcond <- structure(pl$beta, names = colnames(X))
  #FIXME: if we later let glmmTMB.R deal with rank deficient X, then go back to fixef.merMod and copy more complicated part for add.dropped=TRUE case 
  Xzi <- getME(object,"Xzi")
  ffzi <- structure(pl$betazi, names = colnames(Xzi))
  ff=list("conditional model"=ffcond, "zero_inflation"=ffzi)
  l <-sapply(ff, length)>0
  if(sum(l)==1) return(ff[[which(l)]])
  else return(ff[l])
}
##' Extract the modes of the random effects
##'
##' A generic function to extract the conditional modes of the random effects
##' from a fitted model object.  For linear mixed models the conditional modes
##' of the random effects are also the conditional means.
##'
##' If grouping factor i has k levels and j random effects per level the ith
##' component of the list returned by \code{ranef} is a data frame with k rows
##' and j columns.  If \code{condVar} is \code{TRUE} the \code{"postVar"}
##' attribute is an array of dimension j by j by k.  The kth face of this array
##' is a positive definite symmetric j by j matrix.  If there is only one
##' grouping factor in the model the variance-covariance matrix for the entire
##' random effects vector, conditional on the estimates of the model parameters
##' and on the data will be block diagonal and this j by j matrix is the kth
##' diagonal block.  With multiple grouping factors the faces of the
##' \code{"postVar"} attributes are still the diagonal blocks of this
##' conditional variance-covariance matrix but the matrix itself is no longer
##' block diagonal.
##' @name ranef
##' @aliases ranef ranef.glmmTMB
##' @param object an object of a class of fitted models with random effects,
##' typically a \code{"glmmTMB"} object.
##' @param condVar an optional logical argument indicating if the conditional
##' variance-covariance matrices of the random effects should be added as an attribute.
##' @param drop an optional logical argument indicating components of the return
##' value that would be data frames with a single column, usually a column
##' called \sQuote{\code{(Intercept)}}, should be returned as named vectors.
##' @param whichel an optional character vector of names of grouping factors for
##' which the random effects should be returned.  Defaults to all the grouping
##' factors.
##' @param \dots some methods for this generic function require additional
##' arguments.
##' @return A list of data frames, one for each grouping factor for the random
##' effects.  The number of rows in the data frame is the number of levels of
##' the grouping factor.  The number of columns is the dimension of the random
##' effect associated with each level of the factor.
##'
##' If \code{condVar} is \code{TRUE} each of the data frames has an attribute
##' called \code{"postVar"} which is a three-dimensional array with symmetric
##' faces.
##'
##' When \code{drop} is \code{TRUE} any components that would be data frames of
##' a single column are converted to named numeric vectors.
##' @note To produce a \dQuote{caterpillar plot} of the random effects apply
##' \code{\link[lattice:xyplot]{dotplot}} to the result of a call to
##' \code{ranef} with \code{condVar = TRUE}.
##' @examples
##' fm1 <- glmmTMB(Reaction ~ Days + (Days|Subject), sleepstudy)
##' fm2 <- glmmTMB(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
##' fm3 <- glmmTMB(diameter ~ (1|plate) + (1|sample), Penicillin)
##' ranef(fm1)
##' str(rr1 <- ranef(fm1, condVar = TRUE))
##' dotplot(rr1)  ## default
##' ## specify free scales in order to make Day effects more visible
##' dotplot(rr1,scales = list(x = list(relation = 'free')))[["Subject"]]
##' if(FALSE) { ##-- condVar=TRUE is not yet implemented for multiple terms -- FIXME
##' str(ranef(fm2, condVar = TRUE))
##' }
##' op <- options(digits = 4)
##' ranef(fm3, drop = TRUE)
##' options(op)
##' @keywords models methods
##' @importFrom nlme ranef
##' @method ranef glmmTMB
##' @export

ranef.glmmTMB <- function(object, condVar = FALSE, drop = FALSE,
                          whichel = names(ans))
{
  tmp=obj$env$parList()$b
}

##' 
##' Extract or Get Generalize Components from a Fitted Mixed Effects Model
##' Method borrowed from lme4
##' 
##' @importFrom lme4 getME
##' @export getME
##' @method getME glmmTMB
##' @export
getME.glmmTMB <- function(object,
                          name = c("X", "Xzi","Z", "Zzi","Xd","theta"),...)
{
  if(missing(name)) stop("'name' must not be missing")
  stopifnot(is(object,"glmmTMB"))
  ## Deal with multiple names -- "FIXME" is inefficiently redoing things
  if (length(name <- as.character(name)) > 1) {
    names(name) <- name
    return(lapply(name, getME, object = object))
  }
  name <- match.arg(name)
  
  ### Start of the switch
  switch(name,
         "X" = object$obj$env$data$X,
         "Xzi" = object$obj$env$data$Xzi,
         "Z" = object$obj$env$data$Z,
         "Zzi" = object$obj$env$data$Zzi,
         "Xd" = object$obj$env$data$Xd,
         "theta" = object$obj$env$parList()$theta ,
         "..foo.." = # placeholder!
           stop(gettextf("'%s' is not implemented yet",
                         sprintf("getME(*, \"%s\")", name))),
         ## otherwise
         stop(sprintf("Mixed-Effects extraction of '%s' is not available for class \"%s\"",
                      name, class(object))))
}## {getME}

##' Extract the log likelihood of a glmmTMB model
##' 
##' @return object of class \code{logLik} with attributes
##' \item{val}{log likelihood}
##' \item{nobs,nall}{number of non NA observations initially supplied to TMB}
##' \item{df}{number of parameters}
##' @importFrom stats logLik
##' @export logLik
##' @method logLik glmmTMB
##' @export
logLik.glmmTMB<-function(object){
  val <- object$fit$objective
  nobs <- sum(!is.na(object$obj$env$data$yobs))
  structure(val, nobs = nobs, nall = nobs, df = npar.glmmTMB(object),
            class = "logLik")
}

##' Retrieve number of parameters
##' 
##' Also counts dispersion parameter and thetas
npar.glmmTMB<-function(object){
  length(object$fit$par)
}

##' Extracts the variance covariance structure
##' 
##' @importFrom TMB MakeADFun sdreport
##' 
##' @importFrom stats vcov
##' @export vcov
##' @method vcov glmmTMB
##' @export
vcov.glmmTMB<-function(object,full=T){
  if(is.null(object$sdr)){
    warning("Calculating sdreport. Use se=TRUE in glmmTMB to avoid repetitive calculation of sdreport")
    object$sdr<-sdreport(object$obj)
  }

  to_keep <- grep("beta*",colnames(object$sdr$cov.fixed)) # only keep betas
  object$sdr$cov.fixed <- object$sdr$cov.fixed[to_keep,to_keep]
  
  Xnames <- colnames(getME(object,"X"))
  Xzinames <- colnames(getME(object,"Xzi"))
  if(!is.null(Xzinames)) Xzinames <- paste("zi",Xzinames,sep="~")
  Xdnames <- colnames(getME(object,"Xd"))
  if(!is.null(Xdnames)) Xdnames <- paste("d",Xdnames,sep="~")
  colnames(object$sdr$cov.fixed) <- c(Xnames,Xzinames,Xdnames)
  rownames(object$sdr$cov.fixed) <- colnames(object$sdr$cov.fixed)
  
  if(full){
    return(object$sdr$cov.fixed)
  }else{
    di <- grep("d~.*",colnames(object$sdr$cov.fixed))
    zii <- grep("zi~.*",colnames(object$sdr$cov.fixed))
    xi <- setdiff(1:ncol(object$sdr$cov.fixed),c(di,zii))
    output <- list("conditional model" = object$sdr$cov.fixed[xi,xi], "zero_inflation" = object$sdr$cov.fixed[zii,zii], "dispersion" = object$sdr$cov.fixed[di,di])
    l <-sapply(output, length)>0
    if(sum(l)==1) return(output[[which(l)]])
    else return(output[l])
  }
}

