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
fixef.glmmTMB <- function(object, ...) {
  pl <- object$obj$env$parList(object$fit$par, object$obj$env$last.par.best)
  ansc <- structure(pl$beta, names=colnames(getME(object,"X")))
  ansz <- structure(pl$betazi, names=colnames(getME(object,"Xzi")))
  output <- list(conditional_model=ansc, zero_inflation=ansz)
  class(output) <- "fixef.glmmTMB"
  return(output)
}
##' @export
print.fixef.glmmTMB <- function(x, simplify=TRUE, ...) {
  if (simplify && length(x$zero_inflation) == 0L)
    print(unclass(x$conditional_model, ...))
  else
    print(unclass(x), ...)
  invisible(x)
}

##' Extract Random Effects
##'
##' Generic function to extract random effects from \code{glmmTMB} models, both
##' for the conditional model and zero inflation.
##'
##' @param object a \code{glmmTMB} model.
##' @param ... some methods for this generic function require additional
##'   arguments.
##'
##' @return Object of class \code{ranef.glmmTMB} with two components:
##'   \item{conditional_model}{a list of data frames, containing random effects
##'     for the conditional model.}
##'   \item{zero_inflation}{a list of data frames, containing random effects for
##'     the zero inflation.}
##'
##' @note When a model has no zero inflation, the default behavior of
##'   \code{ranef} is to simplify the printed format of the random effects. To
##'   show the full list structure, run \code{print(ranef(model),
##'   simplify=FALSE)}. In all cases, the full list structure is used to access
##'   the data frames (see example).
##'
##' @seealso \code{\link{fixef.glmmTMB}}.
##'
##' @examples
##' data(sleepstudy, package="lme4")
##' model <- glmmTMB(Reaction ~ Days + (1|Subject), sleepstudy)
##' ranef(model)
##' print(ranef(model), simplify=FALSE)
##' ranef(model)$conditional_model$Subject
##'
##' @aliases ranef ranef.glmmTMB
##' @method ranef glmmTMB
##' @export
ranef.glmmTMB <- function(object, ...) {
  ## The arrange() function converts a vector of random effects to a list of
  ## data frames, in the same way as lme4 does.
  arrange <- function(x, listname)
  {
    cnms <- object$modelInfo$reTrms[[listname]]$cnms
    flist <- object$modelInfo$reTrms[[listname]]$flist
    if (!is.null(cnms)) {
      levs <- lapply(fl <- flist, levels)
      asgn <- attr(fl, "assign")
      nc <- vapply(cnms, length, 1L)
      nb <- nc * vapply(levs, length, 1L)[asgn]
      nbseq <- rep.int(seq_along(nb), nb)
      ml <- split(x, nbseq)
      for (i in seq_along(ml))
        ml[[i]] <- matrix(ml[[i]], ncol=nc[i], byrow=TRUE,
                          dimnames=list(NULL, cnms[[i]]))
      x <- lapply(seq_along(fl), function(i)
        data.frame(do.call(cbind, ml[asgn==i]), row.names=levs[[i]],
                   check.names=FALSE))
      names(x) <- names(fl)
    }
    else {
      x <- list()
    }
    return(x)
  }
  pl <- object$obj$env$parList(object$fit$par, object$obj$env$last.par.best)
  ansc <- arrange(pl$b, "condList")
  ansz <- arrange(pl$bzi, "ziList")
  output <- list(conditional_model=ansc, zero_inflation=ansz)
  class(output) <- "ranef.glmmTMB"
  return(output)
}
##' @export
print.ranef.glmmTMB <- function(x, simplify=TRUE, ...) {
  if (simplify && length(x$zero_inflation) == 0L)
    print(unclass(x$conditional_model, ...))
  else
    print(unclass(x), ...)
  invisible(x)
}

##' Extract or Get Generalize Components from a Fitted Mixed Effects Model
##' Method borrowed from lme4
##' 
##' @importFrom lme4 getME
##' @export getME
##' @method getME glmmTMB
##' @export
getME.glmmTMB <- function(object,
                          name = c("X", "Xzi","Z", "Zzi", "Xd", "theta"),
                          ...)
{
  if(missing(name)) stop("'name' must not be missing")
  ## Deal with multiple names -- "FIXME" is inefficiently redoing things
  if (length(name <- as.character(name)) > 1) {
    names(name) <- name
    return(lapply(name, getME, object = object))
  }
  if(name == "ALL") ## recursively get all provided components
      return(sapply(eval(formals()$name),
                    getME.glmmTMB, object=object, simplify=FALSE))

  stopifnot(is(object,"glmmTMB"))
  name <- match.arg(name)
  
  oo.env <- object$obj$env
  ### Start of the switch
  switch(name,
         "X"     = oo.env$data$X,
         "Xzi"   = oo.env$data$Xzi,
         "Z"     = oo.env$data$Z,
         "Zzi"   = oo.env$data$Zzi,
         "Xd"    = oo.env$data$Xd,
         "theta" = oo.env$parList()$theta ,

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
  val <- -object$fit$objective
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
    if(object$modelInfo$allForm$dispformula==~1) object$sdr$cov.fixed <- object$sdr$cov.fixed[-nrow(object$sdr$cov.fixed),-nrow(object$sdr$cov.fixed)]
    
    return(object$sdr$cov.fixed)
  }else{
    di <- grep("d~.*",colnames(object$sdr$cov.fixed))
    zii <- grep("zi~.*",colnames(object$sdr$cov.fixed))
    xi <- setdiff(1:ncol(object$sdr$cov.fixed),c(di,zii))
    output <- list("conditional_model" = object$sdr$cov.fixed[xi,xi], "zero_inflation" = object$sdr$cov.fixed[zii,zii], "dispersion" = object$sdr$cov.fixed[di,di])
    l <-sapply(output, length)>0
    if(sum(l)==1) return(output[[which(l)]])
    else return(output[l])
  }
}

