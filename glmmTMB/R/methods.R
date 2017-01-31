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
  pl <- object$obj$env$parList(object$fit$par, object$fit$parfull)
  structure(list(cond = setNames(pl$beta,   colnames(getME(object, "X"))),
                 zi    = setNames(pl$betazi, colnames(getME(object, "Xzi"))),
                 disp = setNames(pl$betad, colnames(getME(object, "Xd")))),
            class =  "fixef.glmmTMB")
}

## general purpose matching between component names and printable names
cNames <- list(cond = "Conditional model",
               zi = "Zero-inflation model",
               disp = "Dispersion model")

## FIXME: this is a bit ugly. On the other hand, a single-parameter
## dispersion model without a
trivialDisp <- function(object) {
    ## This version works on summary object or fitted model object
    ## FIXME: is there a better way to strip the environment before
    ## comparing?
    identical(deparse(object$call$dispformula),"~1")
}
trivialFixef <- function(xnm,nm) {
    length(xnm)==0 ||
        (nm %in% c('d','disp') && identical(xnm,'(Intercept)'))
    ## FIXME: inconsistent tagging; should change 'Xd' to 'Xdisp'?
}


##' @method print fixef.glmmTMB
##' @export
print.fixef.glmmTMB <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  for(nm in names(x)) {
      if (!trivialFixef(names(x[[nm]]),nm)) {
          cat(sprintf("\n%s:\n", cNames[[nm]]))
          print.default(format(x[[nm]], digits=digits), print.gap = 2L, quote = FALSE)
      }
  }
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

  pl <- getParList(object)
  structure(list(cond = arrange(pl$b, "cond"),
                 zi    = arrange(pl$bzi, "zi")),
            class = "ranef.glmmTMB")
}

##' @method print ranef.glmmTMB
##' @export
print.ranef.glmmTMB <- function(x, simplify=TRUE, ...) {
    print(if (simplify && length(x$zi) == 0L)
              unclass(x$cond) else unclass(x),
          ...)
    invisible(x)
}


##' Extract or Get Generalize Components from a Fitted Mixed Effects Model
##'
##' @aliases getME
##' @param object a fitted \code{glmmTMB} object
##' @param name of the component to be retrieved
##' @param \dots ignored, for method compatibility
##'
##' @seealso \code{\link[lme4]{getME}}
##' Get generic and re-export:
##' @importFrom lme4 getME
##' @export getME
##'
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

## FIXME: (1) why is this non-standard (containing nobs, nall?)
##        (2) do we really need to document it??
## Extract the log likelihood of a glmmTMB model
##
## @return object of class \code{logLik} with attributes
## \item{val}{log likelihood}
## \item{nobs,nall}{number of non NA observations initially supplied to TMB}
## \item{df}{number of parameters}
##' @importFrom stats logLik
##' @export
logLik.glmmTMB <- function(object, ...) {
  if(!is.null(object$sdr)){
    val <- if(object$sdr$pdHess){-object$fit$objective}else{NA}
  }else val <- -object$fit$objective

  nobs <- nobs.glmmTMB(object)
  structure(val, nobs = nobs, nall = nobs, df = length(object$fit$par),
            class = "logLik")
}

##' @importFrom stats nobs
##' @export
nobs.glmmTMB <- function(object, ...) sum(!is.na(object$obj$env$data$yobs))

##' @importFrom stats df.residual
##' @method df.residual glmmTMB
##' @export
##  TODO: not clear whether the residual df should be based
##  on p=length(beta) or p=length(c(theta,beta)) ... but
##  this is just to allow things like aods3::gof to work ...
##  Taken from LME4, including the todo
##
df.residual.glmmTMB <- function(object, ...) {
  nobs(object)-length(object$fit$par)
}


##' Calculate Variance-Covariance Matrix for a Fitted glmmTMB model
##'
##' @param object a \dQuote{glmmTMB} fit
##' @param full return a full variance-covariance matrix?
##' @param \dots ignored, for method compatibility
##' @return By default (\code{full==FALSE}), a list of separate variance-covariance matrices for each model component (conditional, zero-inflation, dispersion).  If \code{full==TRUE}, a single square variance-covariance matrix for \emph{all} top-level model parameters (conditional, dispersion, and variance-covariance parameters)
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats vcov
##' @export
vcov.glmmTMB <- function(object, full=FALSE, ...) {
  if(is.null(sdr <- object$sdr)) {
    warning("Calculating sdreport. Use se=TRUE in glmmTMB to avoid repetitive calculation of sdreport")
    sdr <- sdreport(object$obj)
  }
  keepTag <- if (full) { "."
             } else if (!trivialDisp(object)) { "beta*"
             } else "beta($|[^d])"
  to_keep <- grep(keepTag,colnames(sdr$cov.fixed)) # only keep betas
  covF <- sdr$cov.fixed[to_keep,to_keep,drop=FALSE]

  mkNames <- function(tag) {
      X <- getME(object,paste0("X",tag))
      if (trivialFixef(nn <- colnames(X),tag) &&
          ## if 'full', keep disp even if trivial
          !(full && tag =="d")) character(0)
      else paste(tag,nn,sep="~")
  }

  nameList <- setNames(list(colnames(getME(object,"X")),
                       mkNames("zi"),
                       mkNames("d")),
                names(cNames))
                
  if(full) {
      ## FIXME: haven't really decided if we should drop the
      ##   trivial variance-covariance dispersion parameter ??
      ## if (trivialDisp(object))
      ##    res <- covF[-nrow(covF),-nrow(covF)]

      reNames <- function(tag) {
          re <- object$modelInfo$reStruc[[paste0(tag,"ReStruc")]]
          nn <- mapply(function(n,L) paste(n,seq(L),sep="."),
                 names(re),
                 sapply(re,"[[","blockNumTheta"))
          if (length(nn)==0) return(nn)
          return(paste("theta",gsub(" ","",nn),sep="_"))
      }
      nameList <- c(nameList,list(reNames("cond"),reNames("zi")))

      colnames(covF) <- rownames(covF) <- unlist(nameList)
      res <- covF        ## return just a matrix in this case
  } else {
      splitMat <- function(x) {
          ss <- split(seq_along(colnames(x)),
                      colnames(x))
          lapply(ss,function(z) x[z,z,drop=FALSE])
      }
      covList <- splitMat(covF)
      names(covList) <-
          names(cNames)[match(names(covList),c("beta","betazi","betad"))]
      for (nm in names(covList)) {
          if (length(xnms <- nameList[[nm]])==0) {
              covList[[nm]] <- NULL
          }
          else dimnames(covList[[nm]]) <- list(xnms,xnms)
      }
      res <- covList
      ##  FIXME: should vcov always return a three-element list
      ## (with NULL values for trivial models)?
      class(res) <- c("vcov.glmmTMB","matrix")
  }
  return(res)
}

##' @method print vcov.glmmTMB
##' @export
print.vcov.glmmTMB <- function(x,...) {
    for (nm in names(x)) {
        cat(cNames[[nm]],":\n",sep="")
        print(x[[nm]])
        cat("\n")
    }
    invisible(x)
}

cat.f <- function(...) cat(..., fill = TRUE)

.prt.call.glmmTMB <- function(call, long = TRUE) {
  pass <- 0
  if (!is.null(cc <- call$formula)){
    cat.f("Formula:         ", deparse(cc))
    rhs <- cc[[2]]
    if (!is.null(rhs)) {
        pass<-nchar(deparse(rhs))
    }
  }
  if(!identical(cc <- deparse(call$ziformula),"~0"))
    cat.f("Zero inflation:  ",rep(' ',pass+2), cc, sep='')
  if(!identical(cc <- deparse(call$dispformula),"~1"))
    cat.f("Dispersion:      ",rep(' ',pass+2), cc, sep='')
  if (!is.null(cc <- call$data))
    cat.f("Data:", deparse(cc))
  if (!is.null(cc <- call$weights))
    cat.f("Weights:", deparse(cc))
  if (!is.null(cc <- call$offset))
    cat.f(" Offset:", deparse(cc))
#  if (long && length(cc <- call$control) &&
#      !identical((dc <- deparse(cc)), "lmerControl()"))
    ## && !identical(eval(cc), lmerControl()))
#    cat.f("Control:", dc)
#  if (!is.null(cc <- call$subset))
#    cat.f(" Subset:", deparse(cc))
}


### FIXME: attempted refactoring ...
cat.f2 <- function(call,component,label,lwid,fwid=NULL,cind=NULL) {
    if (!is.null(cc <- call[[component]])) {
        if (!is.null(cind)) {
            ## try to extract component (of formula)
            if (!is.null(ccc <- cc[[cind]]))
                cc <- ccc
        }
        f1 <- format(paste0(label,":"),width=lwid,justify="right")
        f2 <- deparse(cc)
        if (!is.null(fwid)) {
            f2 <- format(f2,width=fwid,justify="right")
        }
        cat(f1,f2,fill=TRUE)
    }
}

## reworked version
.prt.call.glmmTMB2 <- function(call, long = TRUE) {
  labs <- c("Formula","Zero inflation","Dispersion","Data",
            "Weights","Offset","Control","Subset")
  components <- c("formula","ziformula","dispformula",
                  "data","weights","offset","control","subset")

  lwid1 <- max(nchar(labs[1:3]))+2
  for (i in 1:3) {
      cat.f2(call,components[i],labs[i],lwid1,cind=2)
  }
  lwid2 <- max(nchar(labs[-(1:3)]))+1
  for (i in 4:6) {
      cat.f2(call,components[i],labs[i],lwid2)
  }
  if (long && length(cc <- call$control) &&
      (deparse(cc) != "lmerControl()"))
      cat.f2(call,"Control","control",lwid2)
  cat.f2(call,"Subset","subset",lwid2)
}

## following https://github.com/glmmTMB/glmmTMB/issues/134#issuecomment-160805926
## don't use ##' until we're ready to generate a man page
## @param ff name of family (character)
## @param s dispersion (results of sigma(x) for original object
printDispersion <- function(ff,s) {
    ## dispersion
    if (usesDispersion(ff)) {
        if (ff %in% .classicDispersionFamilies) {
            dname <- "Dispersion estimate"
            sname <- "sigma^2"
            sval <- s^2
        } else {
            dname <- "Overdispersion parameter"
            sname <- ""
            sval <- s
        }            
        cat(sprintf("\n%s for %s family (%s): %s",
                    dname,ff,sname,
                    formatC(sval,digits=3)),"\n")
    }
    NULL
}

##' @importFrom lme4 .prt.aictab
##' @method print glmmTMB
##' @export
print.glmmTMB <-
    function(x, digits = max(3, getOption("digits") - 3),
             correlation = NULL, symbolic.cor = FALSE,
             signif.stars = getOption("show.signif.stars"),
             longCall = TRUE, ranef.comp = "Std.Dev.", ...)
{
  ## Type Of Model fit --- REML? ---['class']  & Family & Call
  .prt.call.glmmTMB(x$call, long=longCall)
  ## the 'digits' argument should have an action here
  aictab <- c(AIC = AIC(x), BIC = BIC(x), logLik = logLik(x),
              df.resid = df.residual(x))
  .prt.aictab(aictab, digits=digits+1)
  ## varcorr
  if (!all(sapply(vc <- VarCorr(x),is.null))) {
      cat("Random-effects (co)variances:\n")
      print(VarCorr(x), digits=digits, comp = ranef.comp)
  }
  ## ngroups
  gvec <- list(obs=sprintf("\nNumber of obs: %d",nobs(x)))
  ng <- ngrps.glmmTMB(x)
  for (i in seq_along(ng)) {
      if (length(ng[[i]])>0) {
          nm <- names(ng)[i]
          gvec[[nm]] <- paste0(cNames[nm],": ",
                      paste(paste(names(ng[[i]]), ng[[i]], sep=", "), collapse="; "))
      }
  }
  cat(do.call(paste,c(gvec,list(sep=" / "))),fill=TRUE)

  if(trivialDisp(x)) {# if trivial print here, else below(~x) or none(~0)
    printDispersion(x$modelInfo$familyStr,sigma(x))  
  } 
  ## Fixed effects:
  if(length(cf <- fixef(x)) > 0) {
    cat("\nFixed Effects:\n")
    print(cf, ...)
  } else
    cat("No fixed effect coefficients\n")
  invisible(x)
}

##' @export
model.frame.glmmTMB <- function(formula, ...) {
    formula$frame
}

    
##' Compute residuals for a glmmTMB object
##'
##' @param object a \dQuote{glmmTMB} object
##' @param type (character) residual type
##' @param \dots ignored, for method compatibility
##' @importFrom stats fitted model.response residuals
##' @export
residuals.glmmTMB <- function(object, type=c("response", "pearson"), ...) {
    type <- match.arg(type)
    r <- model.response(object$frame)-fitted(object)
    switch(type,
           response=r,
           pearson={
               if (is.null(v <- family(object)$variance))
                   stop("variance function undefined for family ",
                        sQuote(family(object)$family),"; cannot compute",
                        " Pearson residuals")
               vv <- switch(length(formals(v)),
                            v(fitted(object)),
                            v(fitted(object),sigma(object)),
                            stop("variance function should take 1 or 2 arguments"))
               r/sqrt(vv)
           })
}

## copied from 'stats'

format.perc <- function (probs, digits) {
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), 
    "%")
}

##' @importFrom stats qnorm confint
##' @export
confint.glmmTMB <- function (object, parm, level = 0.95,
                             method=c("Wald","wald",  ## ugh -- allow synonyms?
                                      "profile"),
                             component= "cond", ...) 
{
    dots <- list(...)
    if (length(dots)>0) {
        if (is.null(names(dots))) {
            warning("extra (unnamed) arguments ignored")
        } else {
            warning(paste("extra arguments ignored: ",
                          paste(names(dots),collapse=", ")))
        }
    }
    method <- match.arg(method)
    cf <- unlist(fixef(object)[component])
    pnames <- names(cf)
    if (missing(parm)) 
        parm <- pnames
    else if (is.numeric(parm)) 
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- format.perc(a, 3)
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
        pct))
    if (tolower(method)=="wald") {
        vv <- vcov(object)[component]
        ss <- unlist(lapply(vv,diag))
        ses <- sqrt(ss)[parm]
        ci[] <- cf[parm] + ses %o% fac
    } else {
        stop("profile CI not yet implemented")
        ## FIXME: compute profile(object)
        ## call confint.tmbprofile()
    }
    return(ci)
}

confint.tmbprofile <- function(object, parm=NULL, level = 0.95, ...) {
    ## find locations of top-level (fixed + VarCorr) parameters
    ## fit splines?
    ## invert splines
}

##' @importFrom TMB tmbprofile
profile.glmmTMB <- function(fitted, trace=FALSE, ...) {
    ## lower default spacing?
    ## use Wald std err for initial stepsize guess?
    tmbprofile(fitted$obj, trace=trace, ...)
}

##' @export
## FIXME: establish separate 'terms' components for
##   each model component (conditional, random, zero-inflation, dispersion ...)
terms.glmmTMB <- function(x, component="cond", part="fixed", ...) {
    if (part != "fixed") stop("only fixed terms currently available")
    return(x$modelInfo$reTrms[[component]]$terms[[part]])
    ## terms(x$frame)
}

##' @export
extractAIC.glmmTMB <- function(fit, scale, k = 2, ...) {
    L <- logLik(fit)
    edf <- attr(L,"df")
    return(c(edf,c(-2*L + k*edf)))
}

## deparse(.) returning \bold{one} string
## copied from lme4/R/utilities.R
## Protects against the possibility that results from deparse() will be
##       split after 'width.cutoff' (by default 60, maximally 500)
safeDeparse <- function(x, collapse=" ") paste(deparse(x, 500L), collapse=collapse)

abbrDeparse <- function(x, width=60) {
    r <- deparse(x, width)
    if(length(r) > 1) paste(r[1], "...") else r
}



##' @importFrom methods is
##' @importFrom stats var getCall pchisq anova
##' @export
anova.glmmTMB <- function (object, ..., model.names = NULL) 
{
    mCall <- match.call(expand.dots = TRUE)
    dots <- list(...)
    .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
    ## detect multiple models, i.e. models in ...
    modp <- as.logical(vapply(dots, is, NA, "glmmTMB"))
    if (any(modp)) {
        mods <- c(list(object), dots[modp])
        nobs.vec <- vapply(mods, nobs, 1L)
        if (var(nobs.vec) > 0) 
            stop("models were not all fitted to the same size of dataset")
        if (is.null(mNms <- model.names)) 
            mNms <- vapply(as.list(mCall)[c(FALSE, TRUE, modp)], 
                           safeDeparse, "")
        if (any(duplicated(mNms))) {
            warning("failed to find unique model names, assigning generic names")
            mNms <- paste0("MODEL", seq_along(mNms))
        }
        if (length(mNms) != length(mods)) 
            stop("model names vector and model list have different lengths")
        names(mods) <- sub("@env$", "", mNms)
        llks <- lapply(mods, logLik)
        ii <- order(Df <- vapply(llks, attr, FUN.VALUE = numeric(1), 
            "df"))
        mods <- mods[ii]
        llks <- llks[ii]
        Df <- Df[ii]
        calls <- lapply(mods, getCall)
        data <- lapply(calls, `[[`, "data")
        if (!all(vapply(data, identical, NA, data[[1]]))) 
            stop("all models must be fit to the same data object")
        header <- paste("Data:", abbrDeparse(data[[1]]))
        subset <- lapply(calls, `[[`, "subset")
        if (!all(vapply(subset, identical, NA, subset[[1]]))) 
            stop("all models must use the same subset")
        if (!is.null(subset[[1]])) 
            header <- c(header, paste("Subset:", abbrDeparse(subset[[1]])))
        llk <- unlist(llks)
        chisq <- 2 * pmax(0, c(NA, diff(llk)))
        dfChisq <- c(NA, diff(Df))
        val <- data.frame(Df = Df, AIC = .sapply(llks, AIC), 
            BIC = .sapply(llks, BIC), logLik = llk, deviance = -2 * 
                llk, Chisq = chisq, `Chi Df` = dfChisq, `Pr(>Chisq)` = pchisq(chisq, 
                dfChisq, lower.tail = FALSE), row.names = names(mods), 
            check.names = FALSE)
        class(val) <- c("anova", class(val))
        forms <- lapply(lapply(calls, `[[`, "formula"), deparse)
        structure(val, heading = c(header, "Models:", paste(rep(names(mods), 
            times = lengths(forms)), unlist(forms), sep = ": ")))
    } else stop("no single-model anova() method for glmmTMB")
}

#' @importFrom stats predict
#' @export
fitted.glmmTMB <- function(object, ...) {
    predict(object)
}

.noSimFamilies <- c("beta", "betabinomial", "genpois")

noSim <- function(x) {
    !is.na(match(x, .noSimFamilies))
}

##' Simulate from a glmmTMB fitted model
##' @method simulate glmmTMB 
##' @param object glmmTMB fitted model
##' @param nsim number of response lists to simulate. Defaults to 1.
##' @param seed random number seed
##' @param ... extra arguments 
##' @details Random effects are also simulated from their estimated distribution. 
##' Currently, it is not possible to condition on estimated random effects.  
##' @return returns a list of vectors. The list has length \code{nsim}. 
##' Each simulated vector of observations is the same size as the vector of response variables in the original data set.
##' @importFrom stats simulate
##' @export
simulate.glmmTMB<-function(object, nsim=1, seed=NULL, ...){
    if(noSim(object$modelInfo$family$family))
    {
    	stop("Simulation code has not been implemented for this family")
    }
    if(!is.null(seed)) set.seed(seed)
    ret <- replicate(nsim, object$obj$simulate()$yobs, simplify=FALSE)
    ret
}
