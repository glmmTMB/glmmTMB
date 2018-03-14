##' Extract info from formulas, reTrms, etc., format for TMB
##' @param formula conditional formula
##' @param ziformula zero-inflation formula
##' @param dispformula dispersion formula
##' @param combForm combined formula
##' @param mf call to model frame
##' @param fr model frame
##' @param yobs observed y
##' @param respCol response column
##' @param zioffset offset for zero-inflated model
##' @param doffset offset for dispersion model
##' @param weights weights
##' @param size number of trials in binomial and betabinomial families
##' @param family family object
##' @param se (logical) compute standard error?
##' @param call original \code{glmmTMB} call
##' @param verbose verbosity setting from original \code{glmmTMB} call
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
                       size=NULL,
                       family,
                       se=NULL,
                       call=NULL,
                       verbose=NULL,
                       ziPredictCode="corrected",
                       doPredict=0,
                       whichPredict=integer(0)) {

    ## Handle ~0 dispersion for gaussian family.
    mapArg <- NULL
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
        mapArg <- list(betad = factor(NA)) ## Fix betad
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

    condList  <- getXReTrms(formula, mf, fr)
    ziList    <- getXReTrms(ziformula, mf, fr)
    dispList  <- getXReTrms(dispformula, mf, fr,
                            ranOK=FALSE, "dispersion")

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
    yobs,
    respCol,
    offset = condList$offset,
    zioffset = ziList$offset,
    doffset = dispList$offset,
    weights,
    size,
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
  beta_init <-  if (family$link %in% c("identity","inverse")) 1 else 0

  ## Extra family specific parameters
  numThetaFamily <- (family$family == "tweedie")

  parameters <- with(data.tmb,
                     list(
                       beta    = rep(beta_init, ncol(X)),
                       betazi  = rep(0, ncol(Xzi)),
                       b       = rep(beta_init, ncol(Z)),
                       bzi     = rep(0, ncol(Zzi)),
                       betad   = rep(betad_init, ncol(Xd)),
                       theta   = rep(0, sum(getVal(condReStruc,"blockNumTheta"))),
                       thetazi = rep(0, sum(getVal(ziReStruc,  "blockNumTheta"))),
                       thetaf  = rep(0, numThetaFamily)
                     ))
  randomArg <- c(if(ncol(data.tmb$Z)   > 0) "b",
                 if(ncol(data.tmb$Zzi) > 0) "bzi")
  dispformula <- dispformula.orig ## May have changed - restore
  return(namedList(data.tmb, parameters, mapArg, randomArg, grpVar,
            condList, ziList, dispList, condReStruc, ziReStruc,
            family, respCol,
            allForm=namedList(combForm,formula,ziformula,dispformula),
            fr, se, call, verbose))
}

##' Create X and random effect terms from formula
##' @param formula current formula, containing both fixed & random effects
##' @param mf matched call
##' @param fr full model frame
##' @param ranOK random effects allowed here?
##' @param type label for model type
##' @return a list composed of
##' \item{X}{design matrix for fixed effects}
##' \item{Z}{design matrix for random effects}
##' \item{reTrms}{output from \code{\link{mkReTrms}} from \pkg{lme4}}
##' \item{offset}{offset vector, or vector of zeros if offset not specified}
##'
##' @importFrom stats model.matrix contrasts
##' @importFrom methods new
##' @importFrom lme4 findbars nobars
getXReTrms <- function(formula, mf, fr, ranOK=TRUE, type="") {
    ## fixed-effects model matrix X -
    ## remove random effect parts from formula:
    fixedform <- formula
    RHSForm(fixedform) <- nobars(RHSForm(fixedform))

    nobs <- nrow(fr)
    ## check for empty fixed form

    if (identical(RHSForm(fixedform), ~  0) ||
        identical(RHSForm(fixedform), ~ -1)) {
        X <- NULL
    } else {
        mf$formula <- fixedform
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
        reTrms <- mkReTrms(findbars(RHSForm(formula)), fr)

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
##' blocksize and number of blocks.  Mostly for internal use.
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

##' Fit models with TMB
##' @param formula combined fixed and random effects formula, following lme4
##'     syntax
##' @param data data frame
##' @param family family (variance/link function) information; see \code{\link{family}} for
##' generic family details or \code{\link{family_glmmTMB}} for details of \code{glmmTMB} specific families.  
##' As in \code{\link{glm}}, \code{family} can be specified as (1) a character string
##' referencing an existing family-construction function (e.g. \sQuote{"binomial"}); (2) a symbol referencing
##' such a function (\sQuote{binomial}); or (3) the output of such a function (\sQuote{binomial()}).
##' In addition, for families such as \code{betabinomial} that are special to \code{glmmTMB}, family
##' can be specified as (4) a list comprising the name of the distribution and the link function
##' (\sQuote{list(family="binomial", link="logit")}). However, the first 3 options are preferable.
##' @param ziformula a \emph{one-sided} (i.e., no response variable) formula for
##'     zero-inflation combining fixed and random effects:
##' the default \code{~0} specifies no zero-inflation.
##' Specifying \code{~.} will set the right-hand side of the zero-inflation
##' formula identical to the right-hand side of the main (conditional effects)
##' formula; terms can also be added or subtracted. \strong{Offset terms
##' will automatically be dropped from the conditional effects formula when using \code{~.}}
##' The zero-inflation model uses a logit link.
##' @param dispformula a \emph{one-sided} formula for dispersion containing only fixed effects: the
##'     default \code{~1} specifies the standard dispersion given any family.
##'     The argument is ignored for families that do not have a dispersion parameter.
##'     For an explanation of the dispersion parameter for each family, see (\code{\link{sigma}}).
##'     The dispersion model uses a log link. 
##'     In Gaussian mixed models, \code{dispformula=~0} fixes the parameter to be 0, forcing variance into the random effects.
##' @param weights weights, as in \code{glm}. Not automatically scaled to have sum 1.
##' @param offset offset for conditional model (only):
##' @param se whether to return standard errors
##' @param na.action how to handle missing values (see \code{\link{na.action}})
##' @param verbose logical indicating if some progress indication should be printed to the console.
##' @param doFit whether to fit the full model, or (if FALSE) return the preprocessed data and parameter objects,
##'     without fitting the model
##' @param control control parameters; see \code{\link{glmmTMBControl}}.
##' @importFrom stats gaussian binomial poisson nlminb as.formula terms model.weights
##' @importFrom lme4 subbars findbars mkReTrms nobars
##' @importFrom Matrix t
##' @importFrom TMB MakeADFun sdreport
##' @details
##' \itemize{
##' \item binomial models with more than one trial (i.e., not binary/Bernoulli)
##' can either be specified in the form \code{prob ~ ..., weights = N} or in
##' the more typical two-column matrix (\code{cbind(successes,failures)~...}) form.
##' \item in all cases \code{glmmTMB} returns maximum likelihood estimates - random effects variance-covariance matrices are not REML (so use \code{REML=FALSE} when comparing with \code{lme4::lmer}), and residual standard deviations (\code{\link{sigma}}) are not bias-corrected. Because the \code{\link{df.residual}} method for \code{glmmTMB} currently counts the dispersion parameter, one would need to multiply by \code{sqrt(nobs(fit)/(1+df.residual(fit)))} when comparing with \code{lm} ...
##' \item by default, vector-valued random effects are fitted with
##' unstructured (general positive definite) variance-covariance matrices.
##' Structured variance-covariance matrices can be specified in
##' the form \code{struc(terms|group)}, where \code{struc} is one
##' of
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
##' (note structures marked with * are experimental/untested)
##' }
##' @useDynLib glmmTMB
##' @importFrom stats update
##' @export
##' @examples
##' (m1 <- glmmTMB(count~ mined + (1|site), 
##'   zi=~mined, 
##'   family=poisson, data=Salamanders))
##' summary(m1)
##' \donttest{
##' ## Zero-inflated negative binomial model
##' (m2 <- glmmTMB(count~spp + mined + (1|site), 
##'   zi=~spp + mined, 
##'   family=nbinom2, Salamanders))
##' 
##' ## Hurdle Poisson model
##' (m3 <- glmmTMB(count~spp + mined + (1|site), 
##'   zi=~spp + mined, 
##'   family=truncated_poisson, Salamanders))
##' 
##' ## Binomial model
##' data(cbpp, package="lme4")
##' (tmbm1 <- glmmTMB(cbind(incidence, size-incidence) ~ period + (1 | herd),
##'                data=cbpp, family=binomial))
##' 
##' ## Dispersion model
##' sim1=function(nfac=40, nt=100, facsd=.1, tsd=.15, mu=0, residsd=1)
##' {
##'   dat=expand.grid(fac=factor(letters[1:nfac]), t= 1:nt)
##'   n=nrow(dat)
##'   dat$REfac=rnorm(nfac, sd= facsd)[dat$fac]
##'   dat$REt=rnorm(nt, sd= tsd)[dat$t]
##'   dat$x=rnorm(n, mean=mu, sd=residsd) + dat$REfac + dat$REt
##'   return(dat)
##' }
##' set.seed(101)
##' d1 = sim1(mu=100, residsd =10)
##' d2 = sim1(mu=200, residsd =5)
##' d1$sd="ten"
##' d2$sd="five"
##' dat = rbind(d1, d2)
##' m0 = glmmTMB(x~sd+(1|t), dispformula=~sd, dat)
##' fixef(m0)$disp
##' c(log(5^2), log(10^2)-log(5^2)) #expected dispersion model coefficients
##' }
glmmTMB <- function (
    formula,
    data = NULL,
    family = gaussian(),
    ziformula = ~0,
    dispformula= ~1,
    weights=NULL,
    offset=NULL,
    na.action = getOption("na.action", default=na.fail),
    se=TRUE,
    verbose=FALSE,
    doFit=TRUE,
    control=glmmTMBControl()
    )
{

    ## edited copy-paste from glFormula
    ## glFormula <- function(formula, data=NULL, family = gaussian,
    ##                       subset, weights, na.action, offset,
    ##                       contrasts = NULL, mustart, etastart,
    ##                       control = glmerControl(), ...) {
    call <- mf <- mc <- match.call()

    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
    if (!all(c("family","link") %in% names(family)))
        stop("'family' must contain at least 'family' and 'link' components")
    ## FIXME: warning/message if 'family' doesn't contain 'variance' ?

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
    if (!is.null(offset)) {
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
    etastart <- start <- mustart <- NULL
    if (!is.null(family$initialize)) {
        local(eval(family$initialize))  ## 'local' so it checks but doesn't modify 'y' and 'weights'
    }
    
   if (grepl("^truncated", family$family) &&
       (!is.factor(y) && any(y<1)) & (ziformula == ~0))
        stop(paste0("'", names(respCol), "'", " contains zeros (or values below the allowable range). ",
             "Zeros are compatible with a trucated distribution only when zero-inflation is added."))

    TMBStruc <- 
        mkTMBStruc(formula, ziformula, dispformula,
                   combForm,
                   mf, fr,
                   yobs=y,
                   respCol,
                   weights,
                   family=family,
                   se=se,
                   call=call,
                   verbose=verbose)

    ## Allow for adaptive control parameters
    TMBStruc$control <- lapply(control, eval, envir=TMBStruc)

    ## short-circuit
    if (!doFit) return(TMBStruc)

    ## pack all the bits we will need for fitTMB
    res <- fitTMB(TMBStruc)
    return(res)
}

##' Control parameters for glmmTMB optimization
##' @param optCtrl Passed as argument \code{control} to \code{nlminb}.
##' @param profile Logical; Experimental option to improve speed and
##'                robustness when a model has many fixed effects
##' @param collect Logical; Experimental option to improve speed by
##'                recognizing duplicated observations.
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
##' @export
glmmTMBControl <- function(optCtrl=list(iter.max=300, eval.max=400), 
                           profile=FALSE,
                           collect=FALSE) {
    ## FIXME: Change defaults - add heuristic to decide if 'profile' is beneficial.
    ##        Something like
    ## profile = (length(parameters$beta) >= 2) &&
    ##           (family$family != "tweedie")
    ## (TMB tweedie derivatives currently slow)
    namedList(optCtrl, profile, collect)
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

fitTMB <- function(TMBStruc) {

    control <- TMBStruc$control

    if (control $ collect) {
        ## To avoid side-effects (e.g. nobs.glmmTMB), we restore
        ## original data (with duplicates) after fitting.
        data.tmb.old <- TMBStruc$data.tmb
        TMBStruc$data.tmb <- .collectDuplicates(TMBStruc$data.tmb)
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
        optTime <- system.time(fit <- with(obj,
                                           if( length(par) )
                                               nlminb(start = par, objective = fn, gradient = gr,
                                                      control = control $ optCtrl)
                                           else
                                               list( par=par, objective=fn(par) )
                                           ) )

        sdr <- sdreport(obj, getJointPrecision=TRUE)
        ## FIXME: pdHess can be FALSE
        ##        * Happens for boundary fits (e.g. dispersion close to 0 - see 'spline' example)
        ##        * Option 1: Fall back to old method
        ##        * Option 2: Skip Newton iterations
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
        par <- obj$par; iter <- 0
        if (sdr$pdHess) {
            for (iter in seq_len(5)) { ## FIXME: Make configurable ?
                g <- as.numeric( obj$gr(par) )
                if (max(abs(g)) < 1e-10) break; ## FIXME: Make configurable ?
                par <- par - solve(h, g)
            }
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

    optTime <- system.time(fit <- with(obj, nlminb(start=par, objective=fn,
                                                   gradient=gr, control=control$optCtrl)))
    }

    fit$parfull <- obj$env$last.par.best ## This is in sync with fit$par

    fitted <- NULL

    if (TMBStruc$se) {
        if(control$profile)
            sdr <- sdreport(obj, hessian.fixed=h)
        else
            sdr <- sdreport(obj)
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
                                ## FIXME:apply condList -> cond earlier?
                                reTrms = lapply(list(cond=condList, zi=ziList),
                                                stripReTrms),
                                reStruc = namedList(condReStruc, ziReStruc),
                                allForm))
    ## FIXME: are we including obj and frame or not?
    ##  may want model= argument as in lm() to exclude big stuff from the fit
    ## If we don't include obj we need to get the basic info out
    ##    and provide a way to regenerate it as necessary
    ## If we don't include frame, then we may have difficulty
    ##    with predict() in its current form
    structure(namedList(obj, fit, sdr, call=TMBStruc$call,
                        frame=TMBStruc$fr, modelInfo,
                        fitted),
              class = "glmmTMB")
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
                   family = famL$fami, link = famL$link,
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


