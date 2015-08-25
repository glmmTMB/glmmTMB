##' Extract the fixed-effects estimates
##'
##' Extract the estimates of the fixed-effects parameters from a fitted model.
##' @name fixef
##' @title Extract fixed-effects estimates
##' @aliases fixef fixef.glmmTMB
##' @docType methods
##' @param object any fitted model object from which fixed effects estimates can
##' be extracted.
##' @param \dots optional additional arguments. Currently none are used in any
##' methods.
##' @return a named, numeric vector of fixed-effects estimates.
##' @keywords models
##' @examples
##' data(sleepstudy, package = "lme4")
##' fixef(glmmTMB(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
##' @importFrom nlme fixef
##' @export fixef
##' @export
fixef.glmmTMB <- function(object, ...) {
  pl <- object$obj$env$parList(object$fit$par, object$obj$env$last.par.best)
  structure(list(conditional_model = setNames(pl$beta,   colnames(getME(object, "X"))),
                 zero_inflation    = setNames(pl$betazi, colnames(getME(object, "Xzi")))),
            class =  "fixef.glmmTMB")
}

##' @method print fixef.glmmTMB
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
##' @importFrom nlme ranef
##' @export ranef
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
      x
    }
    else {
      list()
    }
  }

  pl <- object$obj$env$parList(object$fit$par, object$obj$env$last.par.best)
  structure(list(conditional_model = arrange(pl$b, "condList"),
                 zero_inflation    = arrange(pl$bzi, "ziList")),
            class = "ranef.glmmTMB")
}

##' @method print ranef.glmmTMB
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

  stopifnot(inherits(object, "glmmTMB"))
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
##' @export
logLik.glmmTMB <- function(object, ...) {
  val <- -object$fit$objective
  nobs <- sum(!is.na(object$obj$env$data$yobs))
  structure(val, nobs = nobs, nall = nobs, df = npar.glmmTMB(object),
            class = "logLik")
}

##' Retrieve number of parameters
##'
##' Also counts dispersion parameter and thetas
npar.glmmTMB <- function(object){
  length(object$fit$par)
}
## ^^ really ?

##' Extracts the variance covariance structure
##'
##' @importFrom TMB MakeADFun sdreport
##'
##' @importFrom stats vcov
##' @export
vcov.glmmTMB <- function(object, full=TRUE, ...) {
  if(is.null(sdr <- object$sdr)) {
    warning("Calculating sdreport. Use se=TRUE in glmmTMB to avoid repetitive calculation of sdreport")
    sdr <- sdreport(object$obj)
  }

  to_keep <- grep("beta*",colnames(sdr$cov.fixed)) # only keep betas
  covF <- sdr$cov.fixed[to_keep,to_keep]

  Xnames <- colnames(getME(object,"X"))
  Xzinames <- colnames(getME(object,"Xzi"))
  if(!is.null(Xzinames)) Xzinames <- paste("zi",Xzinames,sep="~")
  Xdnames <- colnames(getME(object,"Xd"))
  if(!is.null(Xdnames)) Xdnames <- paste("d",Xdnames,sep="~")
  colnames(covF) <- rownames(covF) <- nms <- c(Xnames,Xzinames,Xdnames)

  if(full) {
    if(object$modelInfo$allForm$dispformula == ~1)
      covF <- covF[-nrow(covF),-nrow(covF)]
    covF
  } else {
    di <- grep("d~.*",nms)
    zii <- grep("zi~.*",nms)
    xi <- setdiff(1:ncol(covF),c(di,zii))
    output <- list(conditional_model = covF[xi,xi],
                   zero_inflation = covF[zii,zii],
                   dispersion = covF[di,di])
    ## FIXME: this may look "convenient" but the result should have same structure in all cases
    l <- sapply(output, length) > 0
    if(sum(l)==1)
      output[[which(l)]]
    else
      output[l]
  }
}

