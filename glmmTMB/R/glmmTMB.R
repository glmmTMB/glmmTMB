##' Extract info from formulas, reTrms, etc., format for TMB
##' @inheritParams glmmTMB
##' @param combForm combined formula
##' @param mf call to model frame
##' @param fr model frame
##' @param yobs observed y
##' @param respCol response column
##' @param zioffset offset for zero-inflated model
##' @param doffset offset for dispersion model
##' @param size number of trials in binomial and betabinomial families
##' @param family family object
##' @param se (logical) compute standard error?
##' @param call original \code{glmmTMB} call
##' @param ziPredictCode zero-inflation code
##' @param doPredict flag to enable sds of predictions
##' @param whichPredict which observations in model frame represent predictions
##' @keywords internal
mkTMBStruc <- function(formula, ziformula, dispformula,
                       combForm,
                       mf, fr,
                       yobs,
                       respCol,
                       ## no conditional offset argument
                       ##  (should be stored in model frame)
                       weights,
                       contrasts,
                       size=NULL,
                       family,
                       se=NULL,
                       call=NULL,
                       verbose=NULL,
                       ziPredictCode="corrected",
                       doPredict=0,
                       whichPredict=integer(0),
                       REML=FALSE,
                       start=NULL,
                       map=NULL) {

  ## handle family specified as naked list
  ## if specified as character or function, should have been converted
  ## to a list with class "family" above ...
    
  ## FIXME: (1) should use proper tryCatch below
  if (!is(family,"family")) {
      ## if family specified as list 
      if (is.list(family)) {
          warning("specifying ",sQuote("family")," as a plain list is deprecated")
      }
      fname <- family$family
      args <- family["link"]
      ff <- try(do.call(fname,args),silent=TRUE)
      if (!inherits(ff,"try-error")) {
          family <- ff
      } else {
          ## fallback: add link information to family
          ## FIXME: is this ever used?
          if (is.null(family$linkfun)) {
              family <- c(family,make.link(family$link))
          }
      }
  }

    ## Handle ~0 dispersion for gaussian family.
    mapArg <- map
    dispformula.orig <- dispformula ## Restore when done
    # family$family contains the *name* of the family
    if ( usesDispersion(family$family) && (dispformula == ~0) ) {
        if (family$family != "gaussian")
            stop("~0 dispersion not implemented for ",
                 sQuote(family$family),
                 " family")
        ## FIXME: Depending on the final estimates, we should somehow
        ## check that this fixed dispersion is small enough.
        betad_init <- log( sqrt( .Machine$double.eps ) )
        dispformula[] <- ~1
        mapArg <- c(mapArg,list(betad = factor(NA))) ## Fix betad
    } else {
        betad_init <- 0
    }

    ## Ignore 'dispformula' argument for non-dispersion families.
    if ( ! usesDispersion(family$family) ) {
        ## if ( dispformula != ~0 &&
        ##      dispformula != ~1 )
        ##     warning(sprintf("dispersion parameter ignored for family %s",
        ##                     sQuote(family)))
        dispformula[] <- ~0
    }

    condList  <- getXReTrms(formula, mf, fr, contrasts=contrasts)
    ziList    <- getXReTrms(ziformula, mf, fr, contrasts=contrasts)
    dispList  <- getXReTrms(dispformula, mf, fr,
                            ranOK=FALSE, type="dispersion",
                            contrasts=contrasts)

    condReStruc <- with(condList, getReStruc(reTrms, ss))
    ziReStruc <- with(ziList, getReStruc(reTrms, ss))

    grpVar <- with(condList, getGrpVar(reTrms$flist))

    nobs <- nrow(fr)

  if (is.null(weights)) weights <- rep(1, nobs)

  ## binomial family:
  ## binomial()$initialize was only executed locally  
  ## yobs could be a factor -> treat as binary following glm
  ## yobs could be cbind(success, failure)
  ## yobs could be binary
  ## (yobs, weights) could be (proportions, size)
  ## On the C++ side 'yobs' must be the number of successes.
  if ( binomialType(family$family) ) {
    if (is.factor(yobs)) {
      ## following glm, ‘success’ is interpreted as the factor not
      ## having the first level (and hence usually of having the
      ## second level).
      yobs <- pmin(as.numeric(yobs)-1,1)
      size <- rep(1, nobs)
    } else {
      if(is.matrix(yobs)) { # yobs=cbind(success, failure)
        size <- yobs[,1] + yobs[,2]
        yobs <- yobs[,1] #successes
      } else {
      if(all(yobs %in% c(0,1))) { #binary
        size <- rep(1, nobs)
      } else { #proportions
          yobs <- weights * yobs
          size <- weights
          weights <- rep(1, nobs)
        }
      }
    }
  }
  if (is.null(size)) size <- numeric(0)

    
  data.tmb <- namedList(
    X = condList$X,
    Z = condList$Z,
    Xzi = ziList$X,
    Zzi = ziList$Z,
    Xd = dispList$X,
    ## Zdisp=dispList$Z,
    ## use c() on yobs, size to strip attributes such as 'AsIs'
    ##  (which confuse MakeADFun)
    yobs = c(yobs),
    respCol,
    offset = condList$offset,
    zioffset = ziList$offset,
    doffset = dispList$offset,
    weights,
    size = c(size),
    ## information about random effects structure
    terms = condReStruc,
    termszi = ziReStruc,
    family = .valid_family[family$family],
    link = .valid_link[family$link],
    ziPredictCode = .valid_zipredictcode[ziPredictCode],
    doPredict = doPredict,
    whichPredict = whichPredict
  )
  getVal <- function(obj, component)
    vapply(obj, function(x) x[[component]], numeric(1))

  ## safer initialization for link functions that might give
  ##  illegal predictions for certain families
  ##  (sqrt() behaves weirdly for beta=0
  ##    [because inverse-link is symmetric around 0?]
  beta_init <-  if (family$link %in% c("identity","inverse","sqrt")) 1 else 0

  ## Extra family specific parameters
  numThetaFamily <- (family$family == "tweedie")

  rr0 <- function(n) {
       if (is.null(n)) numeric(0) else rep(0, n)
  }
  parameters <- with(data.tmb,
                     list(
                       beta    = rep(beta_init, ncol(X)),
                       betazi  = rr0(ncol(Xzi)),
                       b       = rep(beta_init, ncol(Z)),
                       bzi     = rr0(ncol(Zzi)),
                       betad   = rep(betad_init, ncol(Xd)),
                       theta   = rr0(sum(getVal(condReStruc,"blockNumTheta"))),
                       thetazi = rr0(sum(getVal(ziReStruc,  "blockNumTheta"))),
                       thetaf  = rr0(numThetaFamily)
                     ))
    
  for (p in names(start)) {
      if (!(p %in% names(parameters))) {
          stop(sprintf("unrecognized vector '%s' in %s",p,sQuote("start")),
               call. = FALSE)
      }
      if ((Lp <- length(parameters[[p]])) !=  (Ls <- length(start[[p]]))) {
          stop(sprintf("parameter vector length mismatch: in %s, length(%s)==%d, should be %d", sQuote("start"), p, Ls, Lp),
               call. = FALSE)
      }
      parameters[[p]] <- start[[p]]
  }
    
  randomArg <- c(if(ncol(data.tmb$Z)   > 0) "b",
                 if(ncol(data.tmb$Zzi) > 0) "bzi")
  ## REML
  if (REML) randomArg <- c(randomArg, "beta")
  dispformula <- dispformula.orig ## May have changed - restore
  return(namedList(data.tmb, parameters, mapArg, randomArg, grpVar,
            condList, ziList, dispList, condReStruc, ziReStruc,
            family, contrasts, respCol,
            allForm=namedList(combForm,formula,ziformula,dispformula),
            fr, se, call, verbose, REML, map))
}

##' Create X and random effect terms from formula
##' @param formula current formula, containing both fixed & random effects
##' @param mf matched call
##' @param fr full model frame
##' @param ranOK random effects allowed here?
##' @param type label for model type
##' @param contrasts a list of contrasts (see ?glmmTMB)
##' @return a list composed of
##' \item{X}{design matrix for fixed effects}
##' \item{Z}{design matrix for random effects}
##' \item{reTrms}{output from \code{\link{mkReTrms}} from \pkg{lme4}}
##' \item{offset}{offset vector, or vector of zeros if offset not specified}
##'
##' @importFrom stats model.matrix contrasts
##' @importFrom methods new
##' @importFrom lme4 findbars nobars
getXReTrms <- function(formula, mf, fr, ranOK=TRUE, type="", contrasts) {
    ## fixed-effects model matrix X -
    ## remove random effect parts from formula:
    fixedform <- formula
    RHSForm(fixedform) <- nobars(RHSForm(fixedform))

    terms <- NULL ## make sure it's empty in case we don't set it
    
    nobs <- nrow(fr)
    
    ## check for empty fixed form
    ## need to ignore environments when checking!
    ##  ignore.environment= arg only works with closures
    idfun <- function(x,y) {
        environment(x) <- emptyenv()
        environment(y) <- emptyenv()
        return(identical(x,y))
    }
        
    if (idfun(RHSForm(fixedform, as.form=TRUE), ~ 0) ||
        idfun(RHSForm(fixedform, as.form=TRUE), ~ -1)) {
        X <- matrix(ncol=0, nrow=nobs)
        offset <- rep(0,nobs)
    } else {
        tt <- terms(fixedform)
        pv <- attr(mf$formula,"predvars")
        attr(tt, "predvars") <- fix_predvars(pv,tt)
        mf$formula <- tt
        terms_fixed <- terms(eval(mf,envir=environment(fixedform)))
        ## FIXME: make model matrix sparse?? i.e. Matrix:::sparse.model.matrix()
        X <- model.matrix(drop.special2(fixedform), fr, contrasts)
        ## will be 0-column matrix if fixed formula is empty

        offset <- rep(0,nobs)
        terms <- list(fixed=terms(terms_fixed))
        if (inForm(fixedform,quote(offset))) {
            ## hate to match offset terms with model frame names
            ##  via deparse, but since that what was presumably done
            ##  internally to get the model frame names in the first place ...
            for (o in extractForm(fixedform,quote(offset))) {
                offset_nm <- deparse(o)
                ## don't think this will happen, but ...
                if (length(offset_nm)>1) {
                    stop("trouble reconstructing offset name")
                }
                offset <- offset + fr[[offset_nm]]
            }
        }
    }
    
    ## ran-effects model frame (for predvars)
    ## important to COPY formula (and its environment)?
    ranform <- formula
    if (is.null(findbars(ranform))) {
        reTrms <- NULL
        Z <- new("dgCMatrix",Dim=c(as.integer(nobs),0L)) ## matrix(0, ncol=0, nrow=nobs)
        ss <- integer(0)
    } else {

        ## FIXME: check whether predvars are carried along correctly in terms
        if (!ranOK) stop("no random effects allowed in ", type, " term")
        RHSForm(ranform) <- subbars(RHSForm(reOnly(formula)))

        mf$formula <- ranform
        reTrms <- mkReTrms(findbars(RHSForm(formula)), fr, reorder.terms=FALSE)

        ss <- splitForm(formula)
        ss <- unlist(ss$reTrmClasses)

        Z <- t(reTrms$Zt)   ## still sparse ...
    }

    ## if(is.null(rankX.chk <- control[["check.rankX"]]))
    ## rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    ## X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
    ## if(is.null(scaleX.chk <- control[["check.scaleX"]]))
    ##     scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    ## X <- checkScaleX(X, kind=scaleX.chk)

    ## list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula,
    ##      wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))

    namedList(X, Z, reTrms, ss, terms, offset)
}

##' Extract grouping variables for random effect terms from a factor list
##' @title Get Grouping Variable
##' @param x "flist" object; a data frame of factors including an \code{assign} attribute
##' matching columns to random effect terms
##' @return character vector of grouping variables
##' @keywords internal
##' @examples
##' data(cbpp,package="lme4")
##' cbpp$obs <- factor(seq(nrow(cbpp)))
##' rt <- lme4::glFormula(cbind(size,incidence-size)~(1|herd)+(1|obs),
##'   data=cbpp,family=binomial)$reTrms
##' getGrpVar(rt$flist)
##' @export
getGrpVar <- function(x)
{
  assign <- attr(x,"assign")
  names(x)[assign]
}

##' Calculate random effect structure
##' Calculates number of random effects, number of parameters,
##' block size and number of blocks.  Mostly for internal use.
##' @param reTrms random-effects terms list
##' @param ss a character string indicating a valid covariance structure. 
##' Must be one of \code{names(glmmTMB:::.valid_covstruct)};
##' default is to use an unstructured  variance-covariance
##' matrix (\code{"us"}) for all blocks).
##' @return a list
##' \item{blockNumTheta}{number of variance covariance parameters per term}
##' \item{blockSize}{size (dimension) of one block}
##' \item{blockReps}{number of times the blocks are repeated (levels)}
##' \item{covCode}{structure code}
##' @examples
##' data(sleepstudy, package="lme4")
##' rt <- lme4::lFormula(Reaction~Days+(1|Subject)+(0+Days|Subject),
##'                     sleepstudy)$reTrms
##' rt2 <- lme4::lFormula(Reaction~Days+(Days|Subject),
##'                     sleepstudy)$reTrms
##' getReStruc(rt)
##' @importFrom stats setNames dist
##' @export
getReStruc <- function(reTrms, ss=NULL) {

  ## information from ReTrms is contained in cnms, flist
  ## cnms: list of column-name vectors per term
  ## flist: data frame of grouping variables (factors)
  ##   'assign' attribute gives match between RE terms and factors
    if (is.null(reTrms)) {
        list()
    } else {
        ## Get info on sizes of RE components

        assign <- attr(reTrms$flist,"assign")
        nreps <- vapply(assign,
                          function(i) length(levels(reTrms$flist[[i]])),
                          0)
        blksize <- diff(reTrms$Gp) / nreps
        ## figure out number of parameters from block size + structure type

        if (is.null(ss)) {
            ss <- rep("us",length(blksize))
        }

        covCode <- .valid_covstruct[ss]

        parFun <- function(struc, blksize) {
            switch(as.character(struc),
                   "0" = blksize, # diag
                   "1" = blksize * (blksize+1) / 2, # us
                   "2" = blksize + 1, # cs
                   "3" = 2,  # ar1
                   "4" = 2,  # ou
                   "5" = 2,  # exp
                   "6" = 2,  # gau
                   "7" = 3,  # mat
                   "8" = 2 * blksize - 1) # toep
        }
        blockNumTheta <- mapply(parFun, covCode, blksize, SIMPLIFY=FALSE)

        ans <-
            lapply( seq_along(ss), function(i) {
                tmp <-
                    list(blockReps = nreps[i],
                         blockSize = blksize[i],
                         blockNumTheta = blockNumTheta[[i]],
                         blockCode = covCode[i]
                         )
                if(ss[i] == "ar1"){
                    ## FIXME: Keep this warning ?
                    if (any(reTrms$cnms[[i]][1] == "(Intercept)") )
                        warning("AR1 not meaningful with intercept")
                }
                if(ss[i] == "ou"){
                    times <- parseNumLevels(reTrms$cnms[[i]])
                    if (ncol(times) != 1)
                        stop("'ou' structure is for 1D coordinates only.")
                    if (is.unsorted(times, strictly=TRUE))
                        stop("'ou' is for strictly sorted times only.")
                    tmp$times <- drop(times)
                }
                if(ss[i] %in% c("exp", "gau", "mat")){
                    coords <- parseNumLevels(reTrms$cnms[[i]])
                    tmp$dist <- as.matrix( dist(coords) )
                }
                tmp
            })
        setNames(ans, names(reTrms$Ztlist))
    }
}

.noDispersionFamilies <- c("binomial", "poisson", "truncated_poisson")

## BMB: why not just sigma(x)!=1.0 ... ? (redundant with sigma.glmmTMB)
usesDispersion <- function(x) {
    is.na(match(x, .noDispersionFamilies))
    ## !x %in% .noDispersionFamilies
}

.classicDispersionFamilies <- c("gaussian","Gamma","t")

## select only desired pieces from results of getXReTrms
stripReTrms <- function(xrt, whichReTrms = c("cnms","flist"), which="terms") {
  c(xrt$reTrms[whichReTrms],setNames(xrt[which],which))
}

#.okWeightFamilies <- c("binomial", "betabinomial")

okWeights <- function(x) {
	TRUE
  #!is.na(match(x, .okWeightFamilies))
  ## x %in% .okWeightFamilies
}	

## Families for which binomial()$initialize is used
.binomialFamilies <- c("binomial", "betabinomial")
binomialType <- function(x) {
  !is.na(match(x, .binomialFamilies))
}

##' Fit Models with TMB
##'
##' Fit a generalized linear mixed model (GLMM) using Template Model Builder (TMB).
##' @param formula combined fixed and random effects formula, following lme4 syntax.
##' @param data optional data frame containing model variables.
##' @param family a family function, a character string naming a family function, or the result of a call to a family function (variance/link function) information. See \code{\link{family}} for a generic discussion of families or \code{\link{family_glmmTMB}} for details of \code{glmmTMB}-specific families.
##' @param ziformula a \emph{one-sided} (i.e., no response variable) formula for zero-inflation combining fixed and random effects: the default \code{~0} specifies no zero-inflation. Specifying \code{~.} sets the zero-inflation formula identical to the right-hand side of \code{formula} (i.e., the conditional effects formula); terms can also be added or subtracted. \strong{When using \code{~.} as the zero-inflation formula in models where the conditional effects formula contains an offset term, the offset term will automatically be dropped}. The zero-inflation model uses a logit link.
##' @param dispformula a \emph{one-sided} formula for dispersion containing only fixed effects: the default \code{~1} specifies the standard dispersion given any family. The argument is ignored for families that do not have a dispersion parameter. For an explanation of the dispersion parameter for each family, see \code{\link{sigma}}. The dispersion model uses a log link. In Gaussian mixed models, \code{dispformula=~0} fixes the residual variance to be 0 (actually a small non-zero value: at present it is set to \code{sqrt(.Machine$double.eps)}), forcing variance into the random effects.
##' @param weights weights, as in \code{glm}. Not automatically scaled to have sum 1.
##' @param offset offset for conditional model (only).
##' @param contrasts an optional list, e.g., \code{list(fac1="contr.sum")}. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
##' @param na.action how to handle missing values, see \code{\link{na.action}} and \code{\link{model.frame}}. From \code{\link{lm}}: \dQuote{The default is set by the \code{\link{na.action}} setting of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset. The \sQuote{factory-fresh} default is \code{\link{na.omit}}.}
##' @param se whether to return standard errors.
##' @param verbose whether progress indication should be printed to the console.
##' @param doFit whether to fit the full model, or (if FALSE) return the preprocessed data and parameter objects, without fitting the model.
##' @param control control parameters, see \code{\link{glmmTMBControl}}.
##' @param REML whether to use REML estimation rather than maximum likelihood.
##' @param start starting values, expressed as a list with possible components \code{beta}, \code{betazi}, \code{betad} (fixed-effect parameters for conditional, zero-inflation, dispersion models); \code{b}, \code{bzi} (conditional modes for conditional and zero-inflation models); \code{theta}, \code{thetazi} (random-effect parameters, on the standard deviation/Cholesky scale, for conditional and z-i models); \code{thetaf} (extra family parameters, e.g., shape for Tweedie models).
##' @param map a list specifying which parameter values should be fixed to a constant value rather than estimated. \code{map} should be a named list containing factors corresponding to a subset of the internal parameter names (see \code{start} parameter). Distinct factor values are fitted as separate parameter values, \code{NA} values are held fixed: e.g., \code{map=list(beta=factor(c(1,2,3,NA)))} would fit the first three fixed-effect parameters of the conditional model and fix the fourth parameter to its starting value. In general, users will probably want to use \code{start} to specify non-default starting values for fixed parameters. See \code{\link[TMB]{MakeADFun}} for more details.
##' @importFrom stats gaussian binomial poisson nlminb as.formula terms model.weights
##' @importFrom lme4 subbars findbars mkReTrms nobars
##' @importFrom Matrix t
##' @importFrom TMB MakeADFun sdreport
##' @details
##' Binomial models with more than one trial (i.e., not binary/Bernoulli) can either be specified in the form \code{prob ~ ..., weights = N}, or in the more typical two-column matrix \code{cbind(successes,failures)~...} form.
##'
##' Behavior of \code{REML=TRUE} for Gaussian responses matches \code{lme4::lmer}. It may also be useful in some cases with non-Gaussian responses (Millar 2011). Simulations should be done first to verify. 
##'
##' Because the \code{\link{df.residual}} method for \code{glmmTMB} currently counts the dispersion parameter, one would need to multiply by \code{sqrt(nobs(fit) / (1+df.residual(fit)))} when comparing with \code{lm}.
##'
##' By default, vector-valued random effects are fitted with unstructured (general positive definite) variance-covariance matrices. Structured variance-covariance matrices can be specified in the form \code{struc(terms|group)}, where \code{struc} is one of
##' \itemize{
##' \item \code{diag} (diagonal, heterogeneous variance)
##' \item \code{ar1} (autoregressive order-1, homogeneous variance)
##' \item \code{cs} (compound symmetric, heterogeneous variance)
##' \item \code{ou} (* Ornstein-Uhlenbeck, homogeneous variance)
##' \item \code{exp} (* exponential autocorrelation)
##' \item \code{gau} (* Gaussian autocorrelation)
##' \item \code{mat} (* Matérn process correlation)
##' \item \code{toep} (* Toeplitz)
##' }
##' Structures marked with * are experimental/untested.
##'
##' For backward compatibility, the \code{family} argument can also be specified as a list comprising the name of the distribution and the link function (e.g. \code{list(family="binomial", link="logit")}). However, \strong{this alternative is now deprecated}; it produces a warning and will be removed at some point in the future. Furthermore, certain capabilities such as Pearson residuals or predictions on the data scale will only be possible if components such as \code{variance} and \code{linkfun} are present, see \code{\link{family}}.
##'
##' @note
##' For more information about the \pkg{glmmTMB} package, see Brooks et al. (2017) and the \code{vignette(package="glmmTMB")} collection. For the underlying \pkg{TMB} package that performs the model estimation, see Kristensen et al. (2016).
##' @references
##' Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg, C. W., Nielsen, A., Skaug, H. J., Mächler, M. and Bolker, B. M. (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. \emph{The R Journal}, \bold{9}(2), 378--400.
##'
##' Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H. and Bell, B. (2016). TMB: Automatic differentiation and Laplace approximation. \emph{Journal of Statistical Software}, \bold{70}, 1--21.
##'
##' Millar, R. B. (2011). \emph{Maximum Likelihood Estimation and Inference: With Examples in R, SAS and ADMB.} Wiley, New York.
##' @useDynLib glmmTMB
##' @importFrom stats update
##' @export
##' @examples
##' (m1 <- glmmTMB(count ~ mined + (1|site),
##'   zi=~mined,
##'   family=poisson, data=Salamanders))
##' summary(m1)
##' \donttest{
##' ## Zero-inflated negative binomial model
##' (m2 <- glmmTMB(count ~ spp + mined + (1|site),
##'   zi=~spp + mined,
##'   family=nbinom2, data=Salamanders))
##'
##' ## Hurdle Poisson model
##' (m3 <- glmmTMB(count ~ spp + mined + (1|site),
##'   zi=~spp + mined,
##'   family=truncated_poisson, data=Salamanders))
##'
##' ## Binomial model
##' data(cbpp, package="lme4")
##' (bovine <- glmmTMB(cbind(incidence, size-incidence) ~ period + (1|herd),
##'   family=binomial, data=cbpp))
##'
##' ## Dispersion model
##' sim1 <- function(nfac=40, nt=100, facsd=0.1, tsd=0.15, mu=0, residsd=1)
##' {
##'   dat <- expand.grid(fac=factor(letters[1:nfac]), t=1:nt)
##'   n <- nrow(dat)
##'   dat$REfac <- rnorm(nfac, sd=facsd)[dat$fac]
##'   dat$REt <- rnorm(nt, sd=tsd)[dat$t]
##'   dat$x <- rnorm(n, mean=mu, sd=residsd) + dat$REfac + dat$REt
##'   dat
##' }
##' set.seed(101)
##' d1 <- sim1(mu=100, residsd=10)
##' d2 <- sim1(mu=200, residsd=5)
##' d1$sd <- "ten"
##' d2$sd <- "five"
##' dat <- rbind(d1, d2)
##' m0 <- glmmTMB(x ~ sd + (1|t), dispformula=~sd, data=dat)
##' fixef(m0)$disp
##' c(log(5^2), log(10^2)-log(5^2)) # expected dispersion model coefficients
##' }
##'
##' ## Using 'map' to fix random-effects SD to 10
##' m1_map <- update(m1, map=list(theta=factor(NA)),
##'                  start=list(theta=log(10)))
##' VarCorr(m1_map)
glmmTMB <- function(
    formula,
    data = NULL,
    family = gaussian(),
    ziformula = ~0,
    dispformula= ~1,
    weights=NULL,
    offset=NULL,
    contrasts=NULL,
    na.action=na.fail,
    se=TRUE,
    verbose=FALSE,
    doFit=TRUE,
    control=glmmTMBControl(),
    REML=FALSE,
    start=NULL,
    map=NULL
    )
{

    ## edited copy-paste from glFormula
    ## glFormula <- function(formula, data=NULL, family = gaussian,
    ##                       subset, weights, na.action, offset,
    ##                       contrasts = NULL, mustart, etastart,
    ##                       control = glmerControl(), ...) {
    call <- mf <- mc <- match.call()

    if (is.character(family)) {
        if (family=="beta") {
            family <- "beta_family"
            warning("please use ",sQuote("beta_family()")," rather than ",
                    sQuote("\"beta\"")," to specify a Beta-distributed response")
        }
        family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)) {
        ## call family with no arguments
        family <- family()
    }
    ## FIXME: what is this doing? call to a function that's not really
    ##  a family creation function?
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
    fnames <- names(family)
    if (!all(c("family","link") %in% fnames))
        stop("'family' must contain at least 'family' and 'link' components")
    if (length(miss_comp <- setdiff(c("linkfun","variance"),fnames))>0) {
        warning("some components missing from ",sQuote("family"),
                ": downstream methods may fail")
    }
    if (grepl("^quasi", family$family))
        stop('"quasi" families cannot be used in glmmTMB')

    ## extract family and link information from family object
    link <- family$link

    ## lme4 function for warning about unused arguments in ...
    ## ignoreArgs <- c("start","verbose","devFunOnly",
    ##   "optimizer", "control", "nAGQ")
    ## l... <- list(...)
    ## l... <- l...[!names(l...) %in% ignoreArgs]
    ## do.call(checkArgs, c(list("glmer"), l...))

    # substitute evaluated versions
    ## FIXME: denv leftover from lme4, not defined yet

    environment(formula) <- parent.frame()
    call$formula <- mc$formula <- formula
    ## add offset-specified-as-argument to formula as + offset(...)
    ## need evaluate offset within envi
    if (!is.null(eval(substitute(offset),data,
                      enclos=environment(formula)))) {
        formula <- addForm0(formula,makeOp(substitute(offset),op=quote(offset)))
    }


    environment(ziformula) <- environment(formula)
    call$ziformula <- ziformula

    environment(dispformula) <- environment(formula)
    call$dispformula <- dispformula

    ## now work on evaluating model frame
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")

    ## replace . in ziformula with conditional formula, ignoring offset
    if (inForm(ziformula,quote(.))) {
        ziformula <-
            update(RHSForm(drop.special2(formula),as.form=TRUE),
                   ziformula)
    }

    ## want the model frame to contain the union of all variables
    ## used in any of the terms
    ## combine all formulas
    formList <- list(formula, ziformula, dispformula)
    for (i in seq_along(formList)) {
        f <- formList[[i]] ## abbreviate
        ## substitute "|" by "+"; drop specials
        f <- noSpecials(subbars(f),delete=FALSE)
        formList[[i]] <- f
    }
    combForm <- do.call(addForm,formList)
    environment(combForm) <- environment(formula)
    ## model.frame.default looks for these objects in the environment
    ## of the *formula* (see 'extras', which is anything passed in ...),
    ## so they have to be put there ...
    for (i in c("weights", "offset")) {
        if (!eval(bquote(missing(x=.(i)))))
            assign(i, get(i, parent.frame()), environment(combForm))
    }

    mf$formula <- combForm
    fr <- eval(mf,envir=environment(formula),enclos=parent.frame())
    
    ## FIXME: throw an error *or* convert character to factor
    ## convert character vectors to factor (defensive)
    ## fr <- factorize(fr.form, fr, char.only = TRUE)
    ## store full, original formula & offset
    ## attr(fr,"formula") <- combForm  ## unnecessary?
    nobs <- nrow(fr)
    weights <- as.vector(model.weights(fr))

    if(!is.null(weights) & !okWeights(family$family)) {
      stop("'weights' are not available for this family.")
    }

    if (is.null(weights)) weights <- rep(1,nobs)
    
    ## sanity checks (skipped!)
    ## wmsgNlev <- checkNlevels(reTrms$ flist, n=n, control, allow.n=TRUE)
    ## wmsgZdims <- checkZdims(reTrms$Ztlist, n=n, control, allow.n=TRUE)
    ## wmsgZrank <- checkZrank(reTrms$Zt, n=n, control, nonSmall=1e6, allow.n=TRUE)

    ## store info on location of response variable
    respCol <- attr(terms(fr), "response")
    names(respCol) <- names(fr)[respCol]

    ## extract response variable
    ## (name *must* be 'y' to match guts of family()$initialize
    y <- fr[,respCol]
    if (is.matrix(y)) {
        if ( ! binomialType(family$family) ) {
            stop("matrix-valued responses are not allowed")
        }
    }
    
    ## (1) transform 'y' appropriately for binomial models
    ##     (2-column matrix, factor, logical -> numeric)
    ## (2) warn on non-integer values
    ## 'start' should *not* be (reset) to NULL here
    ## as far as we know (i.e. searching src/library/stats/R/family.R
    ##  in the R source), 'start' is only referred to when family="gaussian"
    ##  AND (inverse-link  & any(y==0)) OR (log-link & any(y<=0))
    ##  and then only to check whether it's NULL or not ...
    etastart <- mustart <- NULL

    if (!is.null(family$initialize)) {
        local(eval(family$initialize))  ## 'local' so it checks but doesn't modify 'y' and 'weights'
    }
    
   if (grepl("^truncated", family$family) &&
       (!is.factor(y) && any(y<0.001)) && (ziformula == ~0)) {
        stop(paste0("'", names(respCol), "'", " contains zeros (or values close to zero). ",
             "Zeros are compatible with a truncated distribution only when zero-inflation is added"))
   }

    if (grepl("(nbinom|pois)",family$family)) {
        ## see enum.R: this should cover nbinom1, nbinom2,
        ## poisson, genpois, compois, and the truncated variants
        ## binomial()$initialize already has its own check
        ## (shared by betabinomial)
        if (any(abs(y-round(y))>0.001)) {
            warning(sprintf("non-integer counts in a %s model",
                            family$family))
        }
    }
    
    TMBStruc <- 
        mkTMBStruc(formula, ziformula, dispformula,
                   combForm,
                   mf, fr,
                   yobs=y,
                   respCol,
                   weights,
                   contrasts=contrasts,
                   family=family,
                   se=se,
                   call=call,
                   verbose=verbose,
                   REML=REML,
                   start=start,
                   map=map)

    ## Allow for adaptive control parameters
    TMBStruc$control <- lapply(control, eval, envir=TMBStruc)

    ## short-circuit
    if (!doFit) return(TMBStruc)

    ## pack all the bits we will need for fitTMB
    res <- fitTMB(TMBStruc)
    return(res)
}

##' Control parameters for glmmTMB optimization
##' @param optCtrl   Passed as argument \code{control} to optimizer. Default value (if default \code{nlminb} optimizer is used): \code{list(iter.max=300, eval.max=400)}
##' @param optArgs   additional arguments to be passed to optimizer function (e.g.: \code{list(method="BFGS")} when \code{optimizer=optim})
##' @param profile   Logical; Experimental option to improve speed and
##'                  robustness when a model has many fixed effects
##' @param collect   Logical; Experimental option to improve speed by
##'                  recognizing duplicated observations.
##' @param parallel  Numeric; Setting number of OpenMP threads to evaluate
##'                  the negative log-likelihood in parallel
##' @param optimizer Function to use in model fitting. See \code{Details} for required properties of this function.
##' @importFrom TMB openmp
##' @details
##' The general non-linear optimizer \code{nlminb} is used by
##' \code{\link{glmmTMB}} for parameter estimation. It may sometimes be
##' necessary to tweak some tolerances in order to make a model
##' converge. For instance, the warning \sQuote{iteration limit reached
##' without convergence} may be fixed by increasing the number of
##' iterations using something like
##'
##' \code{glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3))}.
##'
##' The argument \code{profile} allows \code{glmmTMB} to use some special
##' properties of the optimization problem in order to speed up estimation
##' in cases with many fixed effects. Enable this option using
##'
##' \code{glmmTMBControl(profile=TRUE)}.
##'
##' Control parameters may depend on the model specification, because each
##' control component is evaluated inside \code{TMBStruc}, the output
##' of \code{mkTMBStruc}.  To specify that \code{profile} should be
##' enabled for more than 5 fixed effects one can use
##'
##' \code{glmmTMBControl(profile=quote(length(parameters$beta)>=5))}.
##'
##' The \code{optimizer} argument can be any optimization (minimizing) function, provided that:
##' \itemize{
##' \item the first three arguments, in order, are the starting values, objective function, and gradient function;
##' \item it also takes a \code{control} argument;
##' \item it returns a list with elements (at least) \code{convergence} (0 if convergence is successful) and \code{message}
##' }
##' @examples
##' ## fit with default (nlminb) and alternative (optim/BFGS) optimizer
##' m1 <- glmmTMB(count~ mined, family=poisson, data=Salamanders)
##' m1B <- update(m1, control=glmmTMBControl(optimizer=optim,
##'                optArgs=list(method="BFGS")))
##' ## estimates are *nearly* identical:
##' all.equal(fixef(m1), fixef(m1B))
##' @export
glmmTMBControl <- function(optCtrl=NULL,
                           optArgs=list(),
                           optimizer=nlminb,
                           profile=FALSE,
                           collect=FALSE,
                           parallel = NULL) {

    if (is.null(optCtrl) && identical(optimizer,nlminb)) {
        optCtrl <- list(iter.max=300, eval.max=400)
    }
    ## Make sure that we specify at least one thread
    if (!is.null(parallel)) {
        if (is.na(parallel) || parallel < 1) {
            stop("Number of parallel threads must be a numeric >= 1")
        }
        parallel <- as.integer(parallel)
    }
  
    ## FIXME: Change defaults - add heuristic to decide if 'profile' is beneficial.
    ##        Something like
    ## profile = (length(parameters$beta) >= 2) &&
    ##           (family$family != "tweedie")
    ## (TMB tweedie derivatives currently slow)
    namedList(optCtrl, profile, collect, parallel, optimizer, optArgs)
}

##' collapse duplicated observations
##' @keywords internal
##' @importFrom stats runif xtabs
.collectDuplicates <- function(data.tmb) {
    nm <- c("X", "Z", "Xzi", "Zzi", "Xd", "offset",
            "zioffset", "doffset", "yobs",
            "size"[length(data.tmb$size) > 0])
    A <- do.call(cbind, data.tmb[nm])
    ## Restore random seed on exit
    ## FIXME: Simplify ?
    seed <- .GlobalEnv$.Random.seed
    on.exit({
        if (is.null(seed))
            rm(".Random.seed", envir=.GlobalEnv)
        else
            .GlobalEnv$.Random.seed <- seed
    })
    ## Generate hash code for data terms
    hash <- as.vector(A %*% runif(ncol(A)))
    hash <- format(hash, nsmall=20)
    keep <- !duplicated(hash)
    collect <- factor(hash, levels=hash[keep])
    ## Check for collisions
    rownames(A) <- NULL
    A0 <- A[keep, , drop=FALSE]
    if( ! identical (A,
                     A0[unclass(collect), ]) )
        stop("Hash code collision !")
    ## Reduce
    nm <- c("X", "Z", "Xzi", "Zzi", "Xd")
    data.tmb[nm] <- lapply(data.tmb[nm],
                           function(x) x[keep, , drop=FALSE])
    nm <- c("offset", "zioffset", "doffset", "yobs", "size")
    data.tmb[nm] <- lapply(data.tmb[nm],
                           function(x) x[keep])
    ## Update weights
    data.tmb$weights <- xtabs(data.tmb$weights ~ collect)
    data.tmb
}

## FIXME: export fitTMB?
fitTMB <- function(TMBStruc) {

    control <- TMBStruc$control
    
    ## Assign OpenMP threads
    if (!is.null(control$parallel)) {
        n_orig <- TMB::openmp(NULL)
        ## will warn if OpenMP not supported
        ## only proceed farther if OpenMP *is* supported ...
        ## (avoid extra warnings)
        if (n_orig>0) {
            TMB::openmp(n = control$parallel)
            on.exit(TMB::openmp(n = n_orig))
        }
    }

    if (control $ collect) {
        ## To avoid side-effects (e.g. nobs.glmmTMB), we restore
        ## original data (with duplicates) after fitting.
        data.tmb.old <- TMBStruc$data.tmb
        TMBStruc$data.tmb <- .collectDuplicates(TMBStruc$data.tmb)
    }

    ## avoid repetition; rely on environment for parameters
    optfun <- function() {
        with(obj,
             if( length(par) ) {
                 do.call(control$optimizer,
                         c(list(par, fn, gr,
                                control = control $ optCtrl),
                           control $ optArgs))
             } else {
                 list( par=par, objective=fn(par))
             })
    }
    
    if (control $ profile) {
        obj <- with(TMBStruc,
                    MakeADFun(data.tmb,
                              parameters,
                              map = mapArg,
                              random = randomArg,
                              profile = "beta",
                              silent = !verbose,
                              DLL = "glmmTMB"))
        optTime <- system.time(fit <- optfun())

        sdr <- sdreport(obj, getJointPrecision=TRUE)
        parnames <- names(obj$env$par)
        Q <- sdr$jointPrecision; dimnames(Q) <- list(parnames, parnames)
        whichNotRandom <- which( ! parnames %in% c("b", "bzi") )
        Qm <- GMRFmarginal(Q, whichNotRandom)
        h <- as.matrix(Qm) ## Hessian of *all* (non-random) parameters
        TMBStruc$parameters <- obj$env$parList(fit$par, obj$env$last.par.best)
        ## Build object
        obj <- with(TMBStruc,
                    MakeADFun(data.tmb,
                              parameters,
                              map = mapArg,
                              random = randomArg,
                              profile = NULL,
                              silent = !verbose,
                              DLL = "glmmTMB"))
        ## Run up to 5 Newton iterations with fixed (off-mode) hessian
        oldpar <- par <- obj$par; iter <- 0
        ## FIXME: Make configurable ?
        max.newton.steps <- 5
        newton.tol <- 1e-10
        if (sdr$pdHess) {
            ## pdHess can be FALSE
            ##  * Happens for boundary fits (e.g. dispersion close to 0 - see 'spline' example)
            ##    * Option 1: Fall back to old method
            ##    * Option 2: Skip Newton iterations
            for (iter in seq_len(max.newton.steps)) {
                g <- as.numeric( obj$gr(par) )
                if (any(is.na(g)) || max(abs(g)) < newton.tol) break
                par <- par - solve(h, g)
            }
        }
        if (any(is.na(g))) {
            warning("a Newton step failed in profiling")
            par <- oldpar
        }
        fit$par <- par
        fit$objective <- obj$fn(par)
        fit$newton.steps <- iter
    } else {

        obj <- with(TMBStruc,
                    MakeADFun(data.tmb,
                              parameters,
                              map = mapArg,
                              random = randomArg,
                              profile = NULL,
                              silent = !verbose,
                              DLL = "glmmTMB"))
        
        optTime <- system.time(fit <- optfun())
    }

    fit$parfull <- obj$env$last.par.best ## This is in sync with fit$par

    fitted <- NULL

    if (TMBStruc$se) {
        if(control$profile)
            sdr <- sdreport(obj, hessian.fixed=h)
        else
            sdr <- sdreport(obj, getJointPrecision=TMBStruc$REML)
        ## FIXME: assign original rownames to fitted?
    } else {
        sdr <- NULL
    }
    if(!is.null(sdr$pdHess)) {
      if(!sdr$pdHess) {
        warning(paste0("Model convergence problem; ",
                       "non-positive-definite Hessian matrix. ", 
                       "See vignette('troubleshooting')"))
      } else {
        eigval <- try(1/eigen(sdr$cov.fixed)$values, silent=TRUE)
        if( is(eigval, "try-error") || ( min(eigval) < .Machine$double.eps*10 ) ) {
          warning(paste0("Model convergence problem; ",
                       "extreme or very small eigen values detected. ", 
                       "See vignette('troubleshooting')"))
        }
      }
    }

    if ( !is.null(fit$convergence) && fit$convergence != 0)
        warning("Model convergence problem; ",
                fit$message, ". ",
                "See vignette('troubleshooting')")

    if (control $ collect) {
        ## Undo changes made to the data
        TMBStruc$data.tmb <- data.tmb.old
        obj$env$data <- obj$env$dataSanitize(data.tmb.old)
        obj$retape()
    }

    modelInfo <- with(TMBStruc,
                      namedList(nobs=nrow(data.tmb$X),
                                respCol,
                                grpVar,
                                family,
                                contrasts,
                                ## FIXME:apply condList -> cond earlier?
                                reTrms = lapply(list(cond=condList, zi=ziList),
                                                stripReTrms),
                                terms = lapply(list(cond=condList, zi=ziList,
                                                    disp=dispList),
                                               "[[", "terms"),
                                reStruc = namedList(condReStruc, ziReStruc),
                                allForm,
                                REML,
                                map))
    ## FIXME: are we including obj and frame or not?
    ##  may want model= argument as in lm() to exclude big stuff from the fit
    ## If we don't include obj we need to get the basic info out
    ##    and provide a way to regenerate it as necessary
    ## If we don't include frame, then we may have difficulty
    ##    with predict() in its current form

    ret <- structure(namedList(obj, fit, sdr, call=TMBStruc$call,
                        frame=TMBStruc$fr, modelInfo,
                        fitted),
              class = "glmmTMB")

    ## fill in dispersion parameters in environments of family variance
    ## functions, if possible (for glm/effects compatibility)
    ff <- ret$modelInfo$family
    ## family has variance component with extra parameters
    xvarpars <- (length(fv <- ff$variance)>0 &&  
                 length(formals(fv))>1)
    nbfam <- ff$family=="negative.binomial" ||  grepl("nbinom",ff$family)
    if (nbfam || xvarpars) {
        theta <- exp(fit$parfull["betad"]) ## log link
        ## variance() and dev.resids() share an environment
        assign(".Theta",
               theta,
               environment(ret[["modelInfo"]][["family"]][["variance"]]))
    }
    return(ret)
}

##' @importFrom stats AIC BIC
llikAIC <- function(object) {
    llik <- logLik(object)
    AICstats <- 
        c(AIC = AIC(llik), BIC = BIC(llik), logLik = c(llik),
          deviance = -2*llik, ## FIXME:
          df.resid = df.residual(object))
    list(logLik = llik, AICtab = AICstats)
}

## FIXME: export/import from lme4?
ngrps <- function(object, ...) UseMethod("ngrps")

ngrps.default <- function(object, ...) stop("Cannot extract the number of groups from this object")

ngrps.glmmTMB <- function(object, ...) {
    res <- lapply(object$modelInfo$reTrms,
           function(x) vapply(x$flist, nlevels, 1))
    ## FIXME: adjust reTrms names for consistency rather than hacking here
    names(res) <- gsub("List$","",names(res))
    return(res)
    
}

ngrps.factor <- function(object, ...) nlevels(object)


##' @importFrom stats pnorm
##' @method summary glmmTMB
##' @export
summary.glmmTMB <- function(object,...)
{
    if (length(list(...)) > 0) {
        ## FIXME: need testing code
        warning("additional arguments ignored")
    }
    ## figure out useSc
    sig <- sigma(object)

    famL <- family(object)

    mkCoeftab <- function(coefs,vcov) {
        p <- length(coefs)
        coefs <- cbind("Estimate" = coefs,
                       "Std. Error" = sqrt(diag(vcov)))
        if (p > 0) {
            coefs <- cbind(coefs, (cf3 <- coefs[,1]/coefs[,2]),
                           deparse.level = 0)
            ## statType <- if (useSc) "t" else "z"
            statType <- "z"
            ## ??? should we provide Wald p-values???
            coefs <- cbind(coefs, 2*pnorm(abs(cf3), lower.tail = FALSE))
            colnames(coefs)[3:4] <- c(paste(statType, "value"),
                                      paste0("Pr(>|",statType,"|)"))
        }
        coefs
    }

    ff <- fixef(object)
    vv <- vcov(object)
    coefs <- setNames(lapply(names(ff),
            function(nm) if (trivialFixef(names(ff[[nm]]),nm)) NULL else
                             mkCoeftab(ff[[nm]],vv[[nm]])),
                      names(ff))

    llAIC <- llikAIC(object)
                   
    ## FIXME: You can't count on object@re@flist,
    ##	      nor compute VarCorr() unless is(re, "reTrms"):
    varcor <- VarCorr(object)
					# use S3 class for now
    structure(list(logLik = llAIC[["logLik"]],
                   family = famL$family, link = famL$link,
		   ngrps = ngrps(object),
                   nobs = nobs(object),
		   coefficients = coefs, sigma = sig,
		   vcov = vcov(object),
		   varcor = varcor, # and use formatVC(.) for printing.
		   AICtab = llAIC[["AICtab"]], call = object$call
                   ## residuals = residuals(object,"pearson",scaled = TRUE),
		   ## fitMsgs = .merMod.msgs(object),
                   ## optinfo = object@optinfo
		   ), class = "summary.glmmTMB")
               
}

## copied from lme4:::print.summary.merMod (makes use of
##' @importFrom lme4 .prt.family .prt.call .prt.resids .prt.VC .prt.grps
##' @importFrom stats printCoefmat
##' @method print summary.glmmTMB
##' @export
print.summary.glmmTMB <- function(x, digits = max(3, getOption("digits") - 3),
                                 signif.stars = getOption("show.signif.stars"),
                                 ranef.comp = c("Variance", "Std.Dev."),
                                 show.resids = FALSE, ...)
{
    .prt.family(x)
    .prt.call.glmmTMB(x$call); cat("\n")
    .prt.aictab(x$AICtab); cat("\n")
    if (show.resids)
        .prt.resids(x$residuals, digits = digits)

    if (any(whichRE <- !sapply(x$varcor,is.null))) {
        cat("Random effects:\n")
        for (nn in names(x$varcor[whichRE])) {
            cat("\n",cNames[[nn]],":\n",sep="")
            ## lme4:::.prt.VC is not quite what we want here
            print(formatVC(x$varcor[[nn]],
                           digits = digits,
                           comp = ranef.comp),
                  quote=FALSE, digits=digits)
            ## FIXME: redundant nobs output
            .prt.grps(x$ngrps[[nn]],nobs=x$nobs)
        }
    }
    if(trivialDisp(x)) {# if trivial print here, else below(~x) or none(~0)
        printDispersion(x$family,x$sigma)
    }
    for (nn in names(x$coefficients)) {
        cc <- x$coefficients[[nn]]
        p <- length(cc)
        if (p > 0) {
            cat("\n",cNames[[nn]],":\n",sep="")
            printCoefmat(cc, zap.ind = 3, #, tst.ind = 4
                         digits = digits, signif.stars = signif.stars)
        } ## if (p>0)
    }

    invisible(x)
}## print.summary.glmmTMB


