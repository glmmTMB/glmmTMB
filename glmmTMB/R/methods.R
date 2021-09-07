##' Extract Fixed Effects
##'
##' Extract fixed effects from a fitted \code{glmmTMB} model.
##' @name fixef
##' @title Extract fixed-effects estimates
##' @aliases fixef fixef.glmmTMB
##' @docType methods
##' @param object any fitted model object from which fixed effects estimates can
##' be extracted.
##' @param \dots optional additional arguments. Currently none are used in any
##' methods.
##' @return an object of class \code{fixef.glmmTMB} comprising a list of components (\code{cond}, \code{zi}, \code{disp}), each containing a (possibly zero-length) numeric vector of coefficients
##' @keywords models
##' @details The print method for \code{fixef.glmmTMB} object \emph{only displays non-trivial components}: in particular, the dispersion parameter estimate is not printed for models with a single (intercept) dispersion parameter (see examples)
##' @examples
##' data(sleepstudy, package = "lme4")
##' fm1 <- glmmTMB(Reaction ~ Days, sleepstudy)
##' (f1 <- fixef(fm1))
##' f1$cond
##' ## show full coefficients, including dispersion parameter
##' unlist(f1)
##' print.default(f1)
##' @importFrom nlme fixef
##' @export fixef
##' @export
fixef.glmmTMB <- function(object, ...) {
  getXnm <- function(suffix) {
      nm <- paste0("X",suffix)
      return(colnames(getME(object, nm)))
  }
  pl <- object$obj$env$parList(object$fit$par, object$fit$parfull)
  structure(list(cond = setNames(pl$beta,   getXnm("")),
                 zi   = setNames(pl$betazi, getXnm("zi")),
                 disp = setNames(pl$betad,  getXnm("d"))),
            class = "fixef.glmmTMB")
}

## general purpose matching between component names and printable names
cNames <- list(cond = "Conditional model",
               zi = "Zero-inflation model",
               disp = "Dispersion model")

## check identity without worrying about environments etc.
ident <- function(x,target) isTRUE(all.equal(x,target))

formComp <- function(object,type="dispformula",target) {
    ident(object$modelInfo$allForm[[type]],target) ||
        ident(object$call[[type]],target)
}
## FIXME: this is a bit ugly. On the other hand, a single-parameter
## dispersion model without a (... ?)

trivialDisp <- function(object) {
    formComp(object,"dispformula",~1)
}

zeroDisp <- function(object) {
    formComp(object,"dispformula",~0)
}

## no roxygen for now ...
# @param xnm vector of fixed-effect parameter names (e.g. names(fixef(m)$disp))
# @param nm name of component (e.g. "disp")
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
##' Extract random effects from a fitted \code{glmmTMB} model, both
##' for the conditional model and zero inflation.
##'
##' @param object a \code{glmmTMB} model.
##' @param condVar whether to include conditional variances in result.
##' @param \dots some methods for this generic function require additional
##'   arguments (they are unused here and will trigger an error)
##' @return
##' \itemize{
##' \item For \code{ranef}, an object of class \code{ranef.glmmTMB} with two components:
##' \describe{
##'   \item{cond}{a list of data frames, containing random effects
##'     for the conditional model.}
##'   \item{zi}{a list of data frames, containing random effects for
##'     the zero inflation.}
##' }
##' If \code{condVar=TRUE}, the individual list elements within the
##' \code{cond} and \code{zi} components (corresponding to individual
##' random effects terms) will have associated \code{condVar} attributes
##' giving the conditional variances of the random effects values.
##' These are in the form of three-dimensional arrays: see
##' \code{\link{ranef.merMod}} for details. The only difference between
##' the packages is that the attributes are called \sQuote{postVar}
##' in \pkg{lme4}, vs. \sQuote{condVar} in \pkg{glmmTMB}.
##' \item For \code{coef.glmmTMB}: a similar list, but containing
##' the overall coefficient value for each level, i.e., the sum of
##' the fixed effect estimate and the random effect value for that
##' level. \emph{Conditional variances are not yet available as
##' an option for} \code{coef.glmmTMB}.
##' \item For \code{as.data.frame}: a data frame with components
##' \describe{
##' \item{component}{part of the model to which the random effects apply (conditional or zero-inflation)}
##' \item{grpvar}{grouping variable}
##' \item{term}{random-effects term (e.g., intercept or slope)}
##' \item{grp}{group, or level of the grouping variable}
##' \item{condval}{value of the conditional mode}
##' \item{condsd}{conditional standard deviation}
##' }
##' }
##'
##' @note When a model has no zero inflation, the
##' \code{ranef} and \code{coef} print methods simplify the
##' structure shown, by default. To show the full list structure, use
##' \code{print(ranef(model),simplify=FALSE)} or the analogous
##' code for \code{coef}.
##' In all cases, the full list structure is used to access
##' the data frames, see example.
##'
##' @seealso \code{\link{fixef.glmmTMB}}.
##'
##' @examples
##' if (requireNamespace("lme4")) {
##'    data(sleepstudy, package="lme4")
##'    model <- glmmTMB(Reaction ~ Days + (1|Subject), sleepstudy)
##'    rr <- ranef(model)
##'    print(rr, simplify=FALSE)
##'    ## extract Subject conditional modes for conditional model
##'    rr$cond$Subject
##'    as.data.frame(rr)
##' }
##' @aliases ranef ranef.glmmTMB
##' @importFrom nlme ranef
##' @export ranef
##' @export
ranef.glmmTMB <- function(object, condVar=TRUE, ...) {
  check_dots(...)
  ## The arrange() function converts a vector of random effects to a list of
  ## data frames, in the same way as lme4 does.
  ## FIXME: add condVar, make sure format matches lme4
  arrange <- function(x, sd, listname) {
    cnms <- object$modelInfo$reTrms[[listname]]$cnms   ## list of (named) terms and X columns
    reStruc <- object$modelInfo$reStruc[[paste0(listname, "ReStruc")]] ## random-effects structure
    flist <- object$modelInfo$reTrms[[listname]]$flist ## list of grouping variables
    levs <- lapply(flist, levels)
    if (!is.null(cnms)) {  ## FIXME: better test?
      asgn <- attr(flist, "assign")
      ## FIXME: blockReps/blockSize etc. _should_ be stored as integers ...
      nc <- vapply(reStruc, function(x) x$blockSize, numeric(1)) ## number of RE params per block
      nb <- vapply(reStruc, function(x) x$blockReps, numeric(1)) ## number of blocks per RE (may != nlevs in some cases)
      nbseq <- rep.int(seq_along(nb), nb * nc)       ## splitting vector
      ml <- split(x, nbseq)
      for (i in seq_along(ml)) {
          ml[[i]] <- matrix(ml[[i]], ncol = nc[i], byrow = TRUE,
                            dimnames = list(NULL, cnms[[i]]))
      }
      if (!is.null(sd)) {
          sd <- split(sd, nbseq)
          for (i in seq_along(sd)) {
              a <- array(NA, dim=c(nc[i], nc[i], nb[i]))
              ## fill in diagonals: off-diagonals will stay NA (!)
              ## unless we bother to retrieve conditional covariance info
              ## from the fit
              ## when nc>1, what order is the sd vector in?
              ## guessing, level-wise
              for (j in seq(nb[i])) {
                  a[cbind(seq(nc[i]),seq(nc[i]),j)] <-
                      (sd[[i]][nc[i]*(j-1)+seq(nc[i])])^2
              }
              sd[[i]] <- a
          }
      }
      ## combine RE matrices from all terms with the same grouping factor
      x <- lapply(seq_along(flist),
                  function(i) {
                    m <- ml[asgn == i]
                    b2 <- vapply(m, nrow, numeric(1))
                    ub2 <- unique(b2)
                    if (length(ub2)>1)
                      stop("differing numbers of b per group")
                    ## if number of sets of modes != number of levels (e.g. Gaussian process/phyloglmm),
                    ##   generate numeric sequence for names
                    rnms <- if (ub2==length(levs[[i]])) levs[[i]] else seq(ub2)
                    d <- data.frame(do.call(cbind, m),
                               row.names = rnms,
                               check.names = FALSE)

          if (!is.null(sd)) {
              ## attach conditional variance info
              ## called "condVar", *not* "postVar" (contrast to lme4)
              attr(d, "condVar") <- if (length(w <- which(asgn==i))>1) {
                                        ## FIXME: set names?
                                        sd[w]  ## if more than one term, list
                                  } else sd[[w]]  ## else just the array
          }
          return(d)
      })
      names(x) <- names(flist)
      return(x)
    } ## if !is.null(cnms)
    else {
      list()
    }
  } ## arrange()

  pl <- getParList(object)  ## see VarCorr.R
  if (condVar && hasRandom(object))  {
      ss <- summary(object$sdr,"random")
      sdl <- list(b=ss[rownames(ss)=="b","Std. Error"],
                  bzi=ss[rownames(ss)=="bzi","Std. Error"])
  }  else sdl <- NULL
  structure(list(cond = arrange(pl$b, sdl$b, "cond"),
                 zi    = arrange(pl$bzi, sdl$bzi, "zi")),
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

##' @method print coef.glmmTMB
##' @export
print.coef.glmmTMB <- print.ranef.glmmTMB

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
                          name = c("X", "Xzi","Z", "Zzi",
                                   "Xd", "theta", "beta"),
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
  allpars <- oo.env$parList(object$fit$par, object$fit$parfull)
  isSparse <- function(component) { if (is.null(om <- object$modelInfo$sparseX)) FALSE else om[[component]] }
  switch(name,
         "X"     = if (!isSparse("cond")) oo.env$data$X else oo.env$data$XS,
         "Xzi"   = if (!isSparse("zi")) oo.env$data$Xzi else oo.env$data$XziS,
         "Z"     = oo.env$data$Z,
         "Zzi"   = oo.env$data$Zzi,
         "Xd"    = if (!isSparse("disp")) oo.env$data$Xd else oo.env$data$XdS,
         "theta" = allpars$theta ,
         "beta"  = unlist(allpars[c("beta","betazi","betad")]),
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
  df <- sum( ! names(object$fit$parfull) %in% c("b", "bzi") )
  structure(val, nobs = nobs, nall = nobs, df = df,
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
##' @param include_mapped include mapped variables? (these will be given variances and covariances of NA)
##' @param \dots ignored, for method compatibility
##' @return By default (\code{full==FALSE}), a list of separate variance-covariance matrices for each model component (conditional, zero-inflation, dispersion).  If \code{full==TRUE}, a single square variance-covariance matrix for \emph{all} top-level model parameters (conditional, dispersion, and variance-covariance parameters)
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats vcov
##' @export
vcov.glmmTMB <- function(object, full=FALSE, include_mapped=FALSE, ...) {
  check_dots(..., .ignore = "complete")
  REML <- isREML(object)
  if(is.null(sdr <- object$sdr)) {
    warning("Calculating sdreport. Use se=TRUE in glmmTMB to avoid repetitive calculation of sdreport")
    sdr <- sdreport(object$obj, getJointPrecision=REML)
  }
  if (REML) {
      ## NOTE: This code would also work in non-REML case provided
      ## that jointPrecision is present in the object.
      Q <- sdr$jointPrecision
      whichNotRandom <- which( ! rownames(Q) %in% c("b", "bzi") )
      Qm <- GMRFmarginal(Q, whichNotRandom)
      cov.all.parms <- solve(as.matrix(Qm))
  } else {
      cov.all.parms <- sdr$cov.fixed
  }
  keepTag <- if (full) { "."
             } else if (!trivialDisp(object)) { "beta*"
             } else "beta($|[^d])"
  to_keep <- grep(keepTag,colnames(cov.all.parms)) # only keep betas
  covF <- cov.all.parms[to_keep,to_keep,drop=FALSE]

  mkNames <- function(tag) {
      X <- getME(object,paste0("X",tag))
      if (trivialFixef(nn <- colnames(X),tag)
          ## if 'full', keep disp even if trivial, if used by family
          && !(full && tag =="d" &&
               (usesDispersion(family(object)$family) && !zeroDisp(object)))) {
          return(character(0))
      }
      return(paste(tag,nn,sep="~"))
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
      ## nameList for estimated variables;
      nameList <- c(nameList,list(theta=reNames("cond"),thetazi=reNames("zi")))
  }


  ## drop NA-mapped variables

  ## for matching map names vs nameList components ...
  par_components <- c("beta","betazi","betad","theta","thetazi","thetaf")

  fullNameList <- nameList
  map <- object$obj$env$map
  if (length(map)>0) {
        ## fullNameList for all variables, including mapped vars
      ## (nameList will get reduced shortly)
      for (m in seq_along(map)) {
          if (length(NAmap <- which(is.na(map[[m]])))>0) {
              w <- match(names(map)[m],par_components) ##
              if (length(nameList)>=w) { ## may not exist if !full
                  nameList[[w]] <- nameList[[w]][-NAmap]
              }
          }
      }
  }

  if (full) {
        colnames(covF) <- rownames(covF) <- unlist(nameList)
      res <- covF        ## return just a matrix in this case
  } else {
      ## extract block-diagonal matrix
      ss <- split(seq_along(colnames(covF)), colnames(covF))
      covList <- vector("list",3)
      names(covList) <- names(cNames) ## component names
      parnms <- c("beta","betazi", "betad")     ## parameter names
      for (i in seq_along(covList)) {
          nm <- parnms[[i]]
          m <- covF[ss[[nm]],ss[[nm]], drop=FALSE]
          cnm <- names(covList)[[i]]
          xnms <- nameList[[cnm]]
          if (!include_mapped || length(map)==0) {
              dimnames(m) <- list(xnms,xnms)
          } else {
              fnm <- fullNameList[[cnm]]
              mm <- matrix(NA_real_,length(fnm),length(fnm),
                           dimnames=list(fnm,fnm))
              mm[nameList[[cnm]],nameList[[cnm]]] <- m
              m <- mm
          }
          covList[[i]] <- m
      }
      res <- covList[lengths(covList)>0]
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
            dname <- "Dispersion parameter"
            sname <- ""
            sval <- s
        }
        cat(sprintf("\n%s for %s family (%s): %s",
                    dname,ff,sname,
                    formatC(sval,digits=3)),"\n")
    }
    NULL
}

.tweedie_power <- function(object) {
    unname(plogis(get_pars(object)["thetaf"]) + 1)
}

## Print family specific parameters
## @param ff name of family (character)
## @param object glmmTMB output
#' @importFrom stats plogis
printFamily <- function(ff, object) {
    if (ff == "tweedie") {
        power <- .tweedie_power(object)
        cat(sprintf("\nTweedie power parameter: %s",
                    formatC(power, digits=3)), "\n")

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
    printDispersion(x$modelInfo$family$family,sigma(x))
  }
  ## Family specific parameters
  printFamily(x$modelInfo$family$family, x)
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
    if(type=="pearson" &((object$call$ziformula != ~0)|(object$call$dispformula != ~1))) {
        stop("pearson residuals are not implemented for models with zero-inflation or variable dispersion")
    }
    na.act <- attr(object$frame,"na.action")
    mr <- napredict(na.act,model.response(object$frame))
    wts <- model.weights(model.frame(object))
    ## binomial model specified as (success,failure)
    if (!is.null(dim(mr))) {
        wts <- mr[,1]+mr[,2]
        mr <- mr[,1]/wts
    } else if (is.factor(mr)) {
        ## ?binomial:
        ## "‘success’ is interpreted as the factor not having the first level"
        nn <- names(mr)
        mr <- as.numeric(as.numeric(mr)>1)
        names(mr) <- nn  ## restore stripped names
    }
    r <- mr - fitted(object)
    res <- switch(type,
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
               r <- r/sqrt(vv)
               if (!is.null(wts)) {
                   r <- r*sqrt(wts)
               }
               r
    })
    return(res)
}

## Helper to get CI of simple *univariate monotone* parameter
## function, i.e. a function of 'fit$par' and/or 'fit$parfull'.
## Examples: 'sigma.glmmTMB' and some parts of 'VarCorr.glmmTMB'.

## NOT roxygen (##') comments, as we don't want to trigger Rd file creation
## @param object fitted model
## @param f function
## @param reduce
## @param name.prepend
## @param estimate

##' @importFrom stats qchisq
.CI_univariate_monotone <- function(object, f, reduce=NULL,
                                    level=0.95,
                                    name.prepend=NULL,
                                    estimate = TRUE) {
    x <- object
    par <- x$fit$par
    i <- seq_along(x$fit$parfull) ## Pointers into long par vector
    r <- x$obj$env$random
    if(!is.null(r)) i <- i[-r]    ## Pointers into short par subset
    sdr <- x$sdr
    sdpar <- summary(sdr, "fixed")[,2]
    q <- sqrt(qchisq(level, df=1))
    ans <- list()
    x$fit$parfull[i] <- x$fit$par <- par - q * sdpar
    ans$lower <- f(x)
    x$fit$parfull[i] <- x$fit$par <- par + q * sdpar
    ans$upper <- f(x)
    if (estimate) {
        ans$Estimate <- f(object)
    }
    if(is.null(reduce)) reduce <- function(x) x
    ans <- lapply(ans, reduce)
    nm <- names(ans)
    tmp <- cbind(ans$lower, ans$upper)
    if (is.null(tmp) || nrow(tmp) == 0L) return (NULL)
    sort2 <- function(x) if(any(is.na(x))) x*NA else sort(x)
    ans <- cbind( t( apply(tmp, 1, sort2) ) , ans$Estimate )
    colnames(ans) <- nm
    if (!is.null(name.prepend))
        name.prepend <- rep(name.prepend, length.out = nrow(ans))
    rownames(ans) <- paste(name.prepend,
                           rownames(ans), sep="")
    ans
}

## copied from 'stats'

format.perc <- function (probs, digits) {
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
    "%")
}

##' Calculate confidence intervals
##'
##' @details
##' Available methods are
##' \describe{
##' \item{"wald"}{These intervals are based on the standard errors
##' calculated for parameters on the scale
##' of their internal parameterization depending on the family. Derived
##' quantities such as standard deviation parameters and dispersion
##' parameters are back-transformed. It follows that confidence
##' intervals for these derived quantities are typically asymmetric.}
##' \item{"profile"}{This method computes a likelihood profile
##' for the specified parameter(s) using \code{profile.glmmTMB};
##' fits a spline function to each half of the profile; and
##' inverts the function to find the specified confidence interval.}
##' \item{"uniroot"}{This method uses the \code{\link{uniroot}}
##' function to find critical values of one-dimensional profile
##' functions for each specified parameter.}
##' }
##' At present, "wald" returns confidence intervals for variance
##' parameters on the standard deviation/correlation scale, while
##' "profile" and "uniroot" report them on the underlying ("theta")
##' scale: for each random effect, the first set of parameter values
##' are standard deviations on the log scale, while remaining parameters
##' represent correlations on the scaled Cholesky scale (see the
##'
##'
##' @importFrom stats qnorm confint
##' @export
##' @param object \code{glmmTMB} fitted object.
##' @param parm which parameters to profile, specified
#' \itemize{
#' \item by index (position) [\emph{after} component selection for \code{confint}, if any]
#' \item by name (matching the row/column names of \code{vcov(object,full=TRUE)})
#' \item as \code{"theta_"} (random-effects variance-covariance parameters), \code{"beta_"} (conditional and zero-inflation parameters), or \code{"disp_"} or \code{"sigma"} (dispersion parameters)
#' }
#'  Parameter indexing by number may give unusual results when
#'  some parameters have been fixed using the \code{map} argument:
#'  please report surprises to the package maintainers.
##' @param level Confidence level.
##' @param method 'wald', 'profile', or 'uniroot': see Details
##' function)
##' @param component Which of the three components 'cond', 'zi' or
##'     'other' to select. Default is to select 'all'.
##' @param estimate (logical) add a third column with estimate ?
##' @param parallel method (if any) for parallel computation
##' @param ncpus number of CPUs/cores to use for parallel computation
##' @param cl cluster to use for parallel computation
##' @param full CIs for all parameters (including dispersion) ?
##' @param ... arguments may be passed to \code{\link{profile.merMod}} or
##' \code{\link[TMB]{tmbroot}}
##' @examples
##' data(sleepstudy, package="lme4")
##' model <- glmmTMB(Reaction ~ Days + (1|Subject), sleepstudy)
##' model2 <- glmmTMB(Reaction ~ Days + (1|Subject), sleepstudy,
##'     dispformula= ~I(Days>8))
##' confint(model)  ## Wald/delta-method CIs
##' confint(model,parm="theta_")  ## Wald/delta-method CIs
##' confint(model,parm=1,method="profile")
confint.glmmTMB <- function (object, parm = NULL, level = 0.95,
                             method=c("wald",
                                      "Wald",
                                      "profile",
                                      "uniroot"),
                             component = c("all", "cond", "zi", "other"),
                             estimate = TRUE,
                             parallel = c("no", "multicore", "snow"),
                             ncpus = getOption("profile.ncpus", 1L),
                             cl = NULL,
                             full = FALSE,
                             ...) {
    method <- tolower(match.arg(method))
    if (method=="wald") {
        dots <- list(...)
        if (length(dots)>0) {
            if (is.null(names(dots))) {
                warning("extra (unnamed) arguments ignored")
            } else {
                warning(paste("extra arguments ignored: ",
                              paste(names(dots),collapse=", ")))
            }
        }
    }
    components <- match.arg(component, several.ok = TRUE)
    components.has <- function(x)
        any(match(c(x, "all"), components, nomatch=0L)) > 0L
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- format.perc(a, 3)
    fac <- qnorm(a)
    estimate <- as.logical(estimate)
    ci <- matrix(NA, nrow=0, ncol=2 + estimate,
                 dimnames=list(NULL,
                               if (!estimate) pct else c(pct, "Estimate")))

    if (!is.null(parm) || method!="wald") {
        parm <- getParms(parm, object, full)
    }

    wald_comp <- function(component) {
        vv <- vcov(object)[[component]]
        cf <- fixef(object)[[component]]
        ## strip tag (only really necessary for zi~, d~)
        tag <- if (component=="disp") "d" else component
        nn <- gsub(paste0(tag,"~"),"",colnames(vv))
        ## vcov only includes estimated (not mapped/fixed)
        ##  fixed-effect parameters
        cf <- cf[nn]
        ss <- diag(vv)
        ## using [[-extraction; need to add component name explicitly
        ci.tmp <- NULL
        if (length(cf)>0) {
            names(cf) <- names(ss) <-
                paste(component, names(cf), sep=".")
            ses <- sqrt(ss)
            ci.tmp <- cf + ses %o% fac
            if (estimate) ci.tmp <- cbind(ci.tmp, cf)
        }
        return(ci.tmp)
    }

    wald_ci_comp <- function(component) {
        ## VarCorr -> stddev
        cfun <- function(x) {
            ss <- attr(x, "stddev")
            names(ss) <- paste(component,"Std.Dev",names(ss),sep=".")
            cc <- attr(x,"correlation")
            if (length(cc)>1) {
                nn <- outer(colnames(cc),rownames(cc),paste,sep=".")
                cc <- cc[lower.tri(cc)]
                nn <- paste(component,"Cor",nn[lower.tri(nn)],sep=".")
                names(cc) <- nn
                ss <- c(ss,cc)
            }
            return(ss)
        }
        reduce <- function(VC) sapply(VC[[component]], cfun)
        ci.sd <- .CI_univariate_monotone(object,
                                         VarCorr,
                                         reduce = reduce,
                                         level = level,
                                         estimate = estimate)
        ## would consider excluding mapped parameters here
        ## (works automatically for fixed effects via vcov)
        ## but tough because of theta <-> sd/corr mapping;
        ## instead, eliminate rows below where lowerCI==upperCI
        return(ci.sd)
    }

    if (method=="wald") {
        map <- object$modelInfo$map
        for (component in c("cond", "zi") ) {
            if (components.has(component) &&
                length(fixef(object)[[component]])>0) {
                ## variance and estimates
                ci <- rbind(ci, wald_comp(component))
            }
        } ## cond and zi components
        if (components.has("other")) {
            ## sigma
            component <- "disp"
            ff <- object$modelInfo$family$family
            if (usesDispersion(ff)) {
                if (!trivialDisp(object) &&
                    length(fixef(object)[[component]])>0) {
                    ci <- rbind(ci, wald_comp(component))
                } else {
                    ci.sigma <- .CI_univariate_monotone(object,
                                                    sigma,
                                                    reduce = NULL,
                                                    level=level,
                                                    name.prepend="sigma",
                                                    estimate = estimate)
                    ci <- rbind(ci, ci.sigma)
                }
            }
            ## Tweedie power
            if (ff == "tweedie") {
                ci.power <- .CI_univariate_monotone(object,
                                                    .tweedie_power,
                                                    reduce = NULL,
                                                    level=level,
                                                    name.prepend="Tweedie.power",
                                                    estimate = estimate)
                ci <- rbind(ci, ci.power)
            } ## tweedie
        }  ## model has 'other' component
        ## NOW add 'theta' components (match order of params in vcov-full)
        ## FIXME: better to have more robust ordering
        for (component in c("cond", "zi") ) {
            if (components.has(component) &&
                length(ranef(object)[[component]])>0) {
                ci <- rbind(ci, wald_ci_comp(component))
            }
        }

        ## Take subset

        ## drop mapped values (where lower == upper)
        ci <- ci[ci[,2]!=ci[,1], , drop=FALSE]

        ## now get selected parameters
        if (!is.null(parm)) {
            ci <- ci[parm, , drop=FALSE]
        } else {
            ## drop residual std dev/trivial dispersion parameter
            if (!full) {
                ci <- ci[rownames(ci)!="sigma",, drop=FALSE]
            }
        }

        ## end Wald method
    } else if (method=="uniroot") {
        if (isREML(object)) stop("can't compute profiles for REML models at the moment (sorry)")
        ## FIXME: allow greater flexibility in specifying different
        ##  ranges, etc. for different parameters
        plist <- parallel_default(parallel,ncpus)
        parallel <- plist$parallel
        do_parallel <- plist$do_parallel
        FUN <- function(n) {
          n_orig <- openmp(n = object$modelInfo$parallel)
          on.exit(openmp(n_orig))
          TMB::tmbroot(obj=object$obj, name=n, target=0.5*qchisq(level,df=1),
                       ...)
        }
        if (do_parallel) {
            if (parallel == "multicore") {
                L <- parallel::mclapply(parm, FUN, mc.cores = ncpus)
            } else if (parallel=="snow") {
                if (is.null(cl)) {
                    ## start cluster
                    new_cl <- TRUE
                    cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                }
                ## run
                L <- parallel::clusterApply(cl, parm, FUN)
                if (new_cl) {
                    ## stop cluster
                    parallel::stopCluster(cl)
                }
            }
        } else { ## non-parallel
            L <- lapply(as.list(parm), FUN)
        }
        L <- do.call(rbind,L)
        rownames(L) <- rownames(vcov(object,full=TRUE))[parm]
        if (estimate) {
            ee <- object$obj$env
            par <- ee$last.par.best
            if (!is.null(ee$random))
                par <- par[-ee$random]
            par <- par[parm]
            L <- cbind(L,par)
        }
        ci <- rbind(ci,L) ## really just adding column names!
    }
    else {  ## profile CIs
        pp <- profile(object, parm=parm, level_max=level,
                      parallel=parallel,ncpus=ncpus,
                      ...)
        ci <- confint(pp)
    }
    ## if only conditional, strip component prefix
    if (all(substr(rownames(ci),1,5)=="cond.")) {
        rownames(ci) <- sub("^cond\\.","",rownames(ci))
    }
    return(ci)
}

##' @rdname glmmTMB_methods
##' @param x a fitted \code{glmmTMB} object
##' @export
##  modified because e.g. "disp" component didn't get a $reTrms
##  component (updated fitTMB to save a separate "terms" component)
terms.glmmTMB <- function(x, component="cond", part="fixed", ...) {
    if (part != "fixed") stop("only fixed terms currently available")
    if ("terms" %in% names(x$modelInfo)) {
        tt <- x$modelInfo$terms[[component]]
    } else {
        ## allow back-compatibility
        tt <- x$modelInfo$reTrms[[component]]$terms
    }
    return(tt[[part]])
}

##' @export
extractAIC.glmmTMB <- function(fit, scale, k = 2, ...) {
    L <- logLik(fit)
    edf <- attr(L,"df")
    return(c(edf,c(-2*L + k*edf)))
}

## deparse(.) returning \bold{one} string
## previously safeDeparse;
## Protects against the possibility that results from deparse() will be
##       split after 'width.cutoff' (by default 60, maximally 500)
## R >= 4.0.0's deparse1() is a generalization
if((Rv <- getRversion()) < "4.1.0") {
  deparse1 <- function (expr, collapse = " ", width.cutoff = 500L, ...)
      paste(deparse(expr, width.cutoff, ...), collapse = collapse)
}

abbrDeparse <- function(x, width=60) {
    r <- deparse(x, width)
    if(length(r) > 1) paste(r[1], "...") else r
}

sort_termlabs <- function(labs) {
     if (length(labs)==0) return(labs)
     ss <- strsplit(labs, ":")
     ss <- sapply(ss, function(s) { if (length(s)==1) return(s); return(paste(sort(s),collapse=":")) })
     return(sort(ss)) ## sort vector
}

## see whether mod1, mod2 are appropriate for Likelihood ratio testing
CompareFixef <- function (mod1, mod2, component="cond") {
     mr1 <- mod1$modelInfo$REML
     mr2 <- mod2$modelInfo$REML
     if (mr1 != mr2) {
        stop("Can't compare REML and ML fits", call.=FALSE)
     }
     if (mr1 && mr2) {
           tmpf <- function(obj) {   sort_termlabs(attr(terms(obj, component=component),"term.labels")) }
           if (!identical(tmpf(mod1), tmpf(mod2))) {
                stop("Can't compare REML fits with different fixed-effect components", call.=FALSE)
           }
     }
     return(TRUE) ## OK
}

##' @importFrom methods is
##' @importFrom stats var getCall pchisq anova
##' @export
anova.glmmTMB <- function (object, ..., model.names = NULL)
{
    mCall <- match.call(expand.dots = TRUE)
    dots <- list(...)
    ## 'consistent' sapply, i.e. always unlist
    .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
    ## detect multiple models, i.e. models in ...
    modp <- as.logical(vapply(dots, FUN=is, "glmmTMB", FUN.VALUE=NA))
    if (any(modp)) {
        mods <- c(list(object), dots[modp])
        nobs.vec <- vapply(mods, nobs, 1L)
        ## compare all models against first for being fitted consistently;
        ## if all REML, fixed effects must be identical
        vapply(mods[-1], CompareFixef, mod1=mods[[1]], FUN.VALUE=TRUE)
        if (var(nobs.vec) > 0)
            stop("models were not all fitted to the same size of dataset")
        if (is.null(mNms <- model.names))
            mNms <- vapply(as.list(mCall)[c(FALSE, TRUE, modp)],
                           deparse1, "")
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
        ziforms <- lapply(lapply(calls, `[[`, "ziformula"), deparse)
        dispforms <- lapply(lapply(calls, `[[`, "dispformula"), deparse)
        #FIXME only output nontrivial ziforms and dispforms
        structure(val, heading = c(header, "Models:",
            paste(paste(paste(rep(names(mods), times = lengths(forms)), unlist(forms), sep = ": "),
                unlist(ziforms), sep=", zi="),
                unlist(dispforms), sep=", disp=")))
    } else stop("no single-model anova() method for glmmTMB")
}

#' @importFrom stats predict
#' @export
fitted.glmmTMB <- function(object, ...) {
    predict(object,type="response", fast=TRUE)
}

.noSimFamilies <- NULL

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
##' In the binomial family case each simulation is a two-column matrix with success/failure.
##' @importFrom stats simulate
##' @export
simulate.glmmTMB<-function(object, nsim=1, seed=NULL, ...){
    if(noSim(object$modelInfo$family$family))
    {
    	stop("Simulation code has not been implemented for this family")
    }
    ## copied from stats::simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    family <- object$modelInfo$family$family
    ret <- replicate(nsim,
                     object$obj$simulate(par = object$fit$parfull)$yobs,
                     simplify=FALSE)
    if ( binomialType(family) ) {
        size <- object$obj$env$data$size
        ret <- lapply(ret, function(x) cbind(x, size - x, deparse.level=0) )
        class(ret) <- "data.frame"
        rownames(ret) <- as.character(seq_len(nrow(ret[[1]])))
    } else {
        ret <- as.data.frame(ret)
    }
    names(ret) <- paste0("sim_", seq_len(nsim))
    attr(ret, "seed") <- RNGstate
    ret
}

#' Extract the formula of a glmmTMB object
#'
#' @param x a \code{glmmTMB} object
#' @param component formula for which component of the model to return (conditional, zero-inflation, or dispersion)
#' @param fixed.only (logical) drop random effects, returning only the fixed-effect component of the formula?
#' @param ... unused, for generic consistency
#' @importFrom lme4 nobars
#' @export
formula.glmmTMB <- function(x, fixed.only=FALSE,
                            component=c("cond", "zi", "disp"),
                            ...) {
    if (!fixed.only && missing(component)) {
        ## stats::formula.default extracts formula from call
        return(NextMethod(x, ...))
    }
    component <- match.arg(component)
    af <- x$modelInfo$allForm
    ff <- if (component=="cond") af[["formula"]] else af[[paste0(component,"formula")]]
    if (fixed.only) {
        ff <- lme4::nobars(ff)
    }
    return(ff)
}

## need this so we can get contrasts carried through properly to model.matrix


#' Methods for extracting developer-level information from \code{glmmTMB} models
#' @rdname glmmTMB_methods
#' @param object a fitted \code{glmmTMB} object
#' @param component model component ("cond", "zi", or "disp"; not all models contain all components)
#' @param part whether to return results for the fixed or random effect part of the model (at present only \code{part="fixed"} is implemented for most methods)
#' @param \dots additional arguments (ignored or passed to \code{\link{model.frame}})
#' @export

model.matrix.glmmTMB <- function (object, component="cond", part="fixed", ...)
{
    ## FIXME: model.matrix.lm has this stuff -- what does it do/do we want it?
    ## if (n_match <- match("x", names(object), 0L))
    ##    object[[n_match]]
    ## else {
    ## data <- model.frame(object, xlev = object$xlevels, ...)

    ## was calling NextMethod() on the model frame: failed after messing
    ##   with terms structure
    ## could be more efficient to extract $X, $Z rather than re-building
    ## model matrix??
    if (part != "fixed") stop("only fixed model matrices currently available")

    ff <- object$modelInfo$allForm
    form <- ff[[switch(component,
                       cond="formula",
                       zi="ziformula",
                       disp="dispformula")]]
    model.matrix(lme4::nobars(form), model.frame(object, ...),
                 contrasts.arg = object$modelInfo$contrasts)
    ## FIXME: what if contrasts are *different* for different components? (ugh)
    ## should at least write a test to flag this case ...
}

## convert ranef object to a long-format data frame, e.g. suitable
##  for ggplot2 (or homemade lattice plots)
## FIXME: have some gymnastics to do if terms, levels are different
##  for different grouping variables - want to maintain ordering
##  but still allow rbind()ing
##' @export
##' @rdname ranef.glmmTMB
##' @param x a \code{ranef.glmmTMB} object (i.e., the result of running \code{ranef} on a fitted \code{glmmTMB} model)
as.data.frame.ranef.glmmTMB <- function(x, ...) {
    check_dots(..., .ignore = "stringsAsFactors")
    tmpf <- function(x) do.call(rbind,lapply(names(x),asDf0,x=x,id=TRUE))
    x0 <- lapply(x,tmpf)
    x1 <- Map(function(x,n) {
        if (!is.null(x)) x$component <- n; x }, x0, names(x))
    xD <- do.call(rbind,x1)
    ## rename ...
    oldnames <- c("values","ind",".nn","se","id","component")
    newnames <- c("condval","term","grp","condsd","grpvar","component")
    names(xD) <- newnames[match(names(xD),oldnames)]
    ## reorder ...
    neworder <- c("component","grpvar","term","grp","condval")
    if ("condsd" %in% names(xD)) neworder <- c(neworder,"condsd")
    return(xD[neworder])
}

#' @rdname bootmer_methods
#' @title support methods for parametric bootstrapping
#' @param object a fitted glmmTMB object
#' @param newresp a new response vector
#' @export
#' @importFrom lme4 isLMM
#' @importFrom lme4 refit
## don't export refit ...
#' @description see \code{\link[lme4]{refit}} and \code{\link[lme4:isREML]{isLMM}} for details
isLMM.glmmTMB <- function(object) {
   fam <- family(object)
   fam$family=="gaussian" && fam$link=="identity"
}

#' @export
lme4::refit

#' @export
#' @rdname bootmer_methods
#' @importFrom stats formula
#' @param ... additional arguments (for generic consistency; ignored)
#' @examples
#' if (requireNamespace("lme4")) {
#' \dontrun{
#'    fm1 <- glmmTMB(count~mined+(1|spp),
#'                   ziformula=~mined,
#'                   data=Salamanders,
#'                   family=nbinom1)
#'    ## single parametric bootstrap step: refit with data simulated from original model
#'    fm1R <- refit(fm1, simulate(fm1)[[1]])
#'    ## the bootMer function from lme4 provides a wrapper for doing multiple refits
#'    ##   with a specified summary function
#'    b1 <- lme4::bootMer(fm1, FUN=function(x) fixef(x)$zi, nsim=20, .progress="txt")
#'    if (requireNamespace("boot")) {
#'       boot.ci(b1,type="perc")
#'     }
#' }
#' }
#' @details
#' These methods are still somewhat experimental (check your results carefully!), but they should allow parametric bootstrapping.  They work by copying and replacing the original response column in the data frame passed to \code{glmmTMB}, so they will only work properly if (1) the data frame is still available in the environment and (2) the response variable is specified as a single symbol (e.g. \code{proportion} or a two-column matrix constructed on the fly with \code{cbind()}. Untested with binomial models where the response is specified as a factor.
#'
refit.glmmTMB <- function(object, newresp, ...) {
  cc <- getCall(object)
  newdata <- eval.parent(cc$data)
  if (is.null(newdata)) stop("can't locate original 'data' value")
  fresp <- formula(object)[[2]]
  mf0 <- model.frame(object)
  rcol <- attr(attr(mf0, "terms"), "response")
  rnm <- deparse(fresp)
  if (binomialType(family(object)$family)) {
      ## FIXME: check for factor column?
      if ("(weights)" %in% names(mf0)) {
          if (!rnm %in% names(newdata)) stop("can't find response in data")
          w <- rowSums(newresp)
          newdata[[rnm]] <- newresp[,1]/w
          newdata[["(weights)"]] <- w
      } else if (is.matrix(mf0[[rnm]])) {
          if (is.symbol(fresp)) {
              if (!rnm %in% names(newdata)) stop("can't find response in data")
              newdata[[rnm]] <- newresp
          }
          ## matrix response
          else if (identical(quote(cbind),fresp[[1]])) {
              rnm1 <- deparse(fresp[[2]])
              rnm2 <- deparse(fresp[[3]])
              if (!all(c(rnm1,rnm2) %in% names(newdata)))
                  stop("can't find response in data")
              newdata[[rnm1]] <- newresp[,1]
              newdata[[rnm2]] <- newresp[,2]
          } else {
              stop("can't handle this data format, sorry ...")
          }
      } else {
          if (is.matrix(newresp)) newresp <- newresp[,1]
          newdata[[deparse(fresp)]] <- newresp
      }
  } else {
      newdata[[deparse(fresp)]] <- newresp
  }
  cc$data <- quote(newdata)
  return(eval(cc))
}


## copied from lme4, with addition of 'component' argument
## FIXME: migrate back to lme4? component is NULL for back-compat.
## FIXME:
## coef() method for all kinds of "mer", "*merMod", ... objects
## ------  should work with fixef() + ranef()  alone
coefMer <- function(object, component=NULL, ...)
{
    if (length(list(...)))
        warning('arguments named "', paste(names(list(...)), collapse = ", "),
                '" ignored')
    fef <- fixef(object)
    if (!is.null(component)) fef <- fef[[component]]
    fef <- data.frame(rbind(fef), check.names = FALSE)
    ref <- ranef(object)
    if (!is.null(component)) ref <- ref[[component]]
    ## check for variables in RE but missing from FE, fill in zeros in FE accordingly
    refnames <- unlist(lapply(ref,colnames))
    nmiss <- length(missnames <- setdiff(refnames,names(fef)))
    if (nmiss > 0) {
        fillvars <- setNames(data.frame(rbind(rep(0,nmiss))),missnames)
        fef <- cbind(fillvars,fef)
    }
    val <- lapply(ref, function(x)
                  fef[rep.int(1L, nrow(x)),,drop = FALSE])
    for (i in seq(a = val)) {
        refi <- ref[[i]]
        row.names(val[[i]]) <- row.names(refi)
        nmsi <- colnames(refi)
        if (!all(nmsi %in% names(fef)))
            stop("unable to align random and fixed effects")
        for (nm in nmsi) val[[i]][[nm]] <- val[[i]][[nm]] + refi[,nm]
    }
    class(val) <- "coef.mer"
    val
} ##  {coefMer}

#' @rdname ranef.glmmTMB
#' @export
coef.glmmTMB <- function(object,
                         condVar=FALSE, ...) {
    model.has.component <- function(x) {
        !is.null(object$modelInfo$reTrms[[x]]$cnms)
    }
    get.coef <- function(x) {
        if (!model.has.component(x)) return(list())
        return(coefMer(object, component=x))
    }
    res <- list(
        cond = get.coef("cond"),
        zi = get.coef("zi")
    )
    if (condVar) {
        stop("condVar not (yet) available for coefficients")
        sdr <- TMB::sdreport(object$obj, getJointPrecision=TRUE)
        v <- solve(sdr$jointPrecision)
        ## FIXME:: sort out variance calculation, using Z and X
    }
    class(res) <- "coef.glmmTMB"
    return(res)
}


##' Extract weights from a glmmTMB object
##'
##' @details
##' At present only explicitly specified
##' \emph{prior weights} (i.e., weights specified
##' in the \code{weights} argument) can be extracted from a fitted model.
##' \itemize{
##' \item Unlike other GLM-type models such as \code{\link{glm}} or
##' \code{\link[lme4]{glmer}}, \code{weights()} does not currently return
##' the total number of trials when binomial responses are specified
##' as a two-column matrix.
##' \item Since \code{glmmTMB} does not fit models via iteratively
##' weighted least squares, \code{working weights} (see \code{\link[stats:glm]{weights.glm}}) are unavailable.
##' }
##' @importFrom stats model.frame
##' @importFrom stats weights
##' @param object a fitted \code{glmmTMB} object
##' @param type weights type
##' @param ... additional arguments (not used; for methods compatibility)
##' @export
weights.glmmTMB <- function(object, type="prior", ...) {
    type <- match.arg(type)  ## other types are *not* OK
    if (length(list(...)>0)) {
        warning("unused arguments ignored: ",
             paste(shQuote(names(list(...))),collapse=","))
    }
    stats::model.frame(object)[["(weights)"]]
}

# would like to export this only as a method, but not sure how ...
# https://stackoverflow.com/questions/29079179/does-using-package-generics-require-the-package-to-be-in-depends-or-imports

# extract model parameters
#
# This is a utility function for multcomp::glht
#
## @param model fitted glmmTMB model
## @param coef. function for retrieving coefficients
## @param vcov. function for retrieving covariance matrix
## @param df degrees of freedom
## @param component which model component to test (cond, zi, or disp)

##' @rawNamespace if(getRversion() >= "3.6.0") {
##'      S3method(multcomp::modelparm, glmmTMB)
##' } else {
##'    export(modelparm.glmmTMB)
##' }
modelparm.glmmTMB <- function (model, coef. = function(x) fixef(x)[[component]],
                               vcov. = function(x) vcov(x)[[component]],
                               df = NULL, component="cond", ...) {
    beta <- coef.(model)
    sigma <- vcov.(model)
    estimable <- unname(!is.na(beta))
    if (is.null(df)) df <- 0
    RET <- list(coef = beta, vcov = sigma, df = df, estimable = estimable)
    class(RET) <- "modelparm"
    return(RET)
}
