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
  structure(list(cond = setNames(pl$beta,   colnames(getME(object, "X"))),
                 zi    = setNames(pl$betazi, colnames(getME(object, "Xzi"))),
                 disp = setNames(pl$betad, colnames(getME(object, "Xd")))),
            class =  "fixef.glmmTMB")
}

##' @method print fixef.glmmTMB
##' @export
print.fixef.glmmTMB<-function(object, digits = max(3, getOption("digits") - 3)){
  name <- list(cond = "Conditional model", zi = "Zero inflation", disp = "Dispersion")
  for(x in names(object)){
    if((length(object[[x]])-as.numeric(x == 'disp' & '(Intercept)' %in% names(object[[x]])))>0){
      cat("\n",name[[x]],"\n",sep="")
      print.default(format(object[[x]],digits=digits), print.gap = 2L, quote = FALSE)
    }
  }
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
  structure(list(cond = arrange(pl$b, "condList"),
                 zi    = arrange(pl$bzi, "ziList")),
            class = "ranef.glmmTMB")
}

##' @method print ranef.glmmTMB
##' @export
print.ranef.glmmTMB <- function(x, simplify=TRUE, ...) {
  if (simplify && length(x$zi) == 0L)
    print(unclass(x$cond, ...))
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
  structure(val, nobs = nobs, nall = nobs, df = length(object$fit$par),
            class = "logLik")
}

##' @importFrom stats nobs
##' @method nobs glmmTMB
nobs.glmmTMB <- function(object, ...) sum(!is.na(object$obj$env$data$yobs))


.prt.aictab <- function(object, digits = 1) {
  aictab <- c(AIC = AIC(object), BIC = BIC(object), logLik = logLik(object),
              df.resid = df.residual(object))
  t.4 <- round(aictab, digits)

    ## slight hack to get residual df formatted as an integer
    t.4F <- format(t.4)
    t.4F["df.resid"] <- format(t.4["df.resid"])
    print(t.4F, quote = FALSE)

}

##' Residual Degrees-of-Freedom
##'
##' @importFrom stats df.residual
##' @method df.residual glmmTMB
##  TODO: not clear whether the residual df should be based
##  on p=length(beta) or p=length(c(theta,beta)) ... but
##  this is just to allow things like aods3::gof to work ...
##  Taken from LME4, including the todo
##
df.residual.glmmTMB <- function(object, ...) {
  nobs(object)-length(object$fit$par)
}


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

cat.f <- function(...) cat(..., fill = TRUE)

.prt.call.glmmTMB <- function(call, long = TRUE) {
  pass<-0
  if (!is.null(cc <- call$formula)){
    cat.f("Formula:         ", deparse(cc))
    pass<-nchar(as.character(call$formula[[2]]))
  }
  if(!is.null(cc <- call$ziformula))
    cat.f("Zero inflation:  ",rep(' ',pass+2),'~ ' ,deparse(cc[[2]]),sep='')
  if(!is.null(cc <- call$dispformula))
    cat.f("Dispersion:      ",rep(' ',pass+2),'~ ', deparse(cc[[2]]), sep='')
  if (!is.null(cc <- call$data))
    cat.f("   Data:", deparse(cc))
  if (!is.null(cc <- call$weights))
    cat.f("Weights:", deparse(cc))
  if (!is.null(cc <- call$offset))
    cat.f(" Offset:", deparse(cc))
  if (long && length(cc <- call$control) &&
      !identical((dc <- deparse(cc)), "lmerControl()"))
    ## && !identical(eval(cc), lmerControl()))
    cat.f("Control:", dc)
  if (!is.null(cc <- call$subset))
    cat.f(" Subset:", deparse(cc))
}

##' Print glmmTMB model
##' @method print glmmTMB
##' @export
##' 
print.glmmTMB<-function(object, digits = max(3, getOption("digits") - 3),
                        correlation = NULL, symbolic.cor = FALSE,
                        signif.stars = getOption("show.signif.stars"),
                        ranef.comp = "Std.Dev.", ...){

  # TYPE OF MODEL FIT --- REML? ---['class']
  # FAMILY
  # CALL
  .prt.call.glmmTMB(object$call)
  # AIC TABLE
  .prt.aictab(object,4)
  # varcorr
  # ngroups
  
  # Print fixed effects
  if(length(cf <- fixef(object)) > 0) {
    cat("\nFixed Effects:\n")
    print(cf, ...)
  } else cat("No fixed effect coefficients\n")
  
}
