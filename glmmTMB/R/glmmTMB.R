## internal flag for debugging OpenMP behaviour
debug_openmp <- FALSE

## glmmTMB openmp controller copied from TMB (Windows needs it).
openmp <- function (n = NULL) {
    if (debug_openmp && !is.null(n)) {
        cat("setting OpenMP threads to ", n, "\n")
    }
    ## FIXME: redundant with integer-setting within omp_num_threads C++ def in utils.cpp
    null_arg <- is.null(n)
    if (!null_arg) n <- as.integer(n)
    ## only want to warn if attempt to set >1 threads in absence
    ## of OpenMP support ..
    if (null_arg || n <= 1) {
      w <- options(warn = -1)
      on.exit(options(warn = w[["warn"]]))
    }
    TMB::openmp(n, DLL="glmmTMB")
}

##' Change starting parameters, either by residual method or by user input (start)
##' @inheritParams mkTMBStruc
##' @param formula current formula, containing both fixed & random effects
##' @param ziformula a \emph{one-sided} (i.e., no response variable) formula for zero-inflation combining fixed and random effects: the default \code{~0} specifies no zero-inflation. Specifying \code{~.} sets the zero-inflation formula identical to the right-hand side of \code{formula} (i.e., the conditional effects formula); terms can also be added or subtracted. \strong{When using \code{~.} as the zero-inflation formula in models where the conditional effects formula contains an offset term, the offset term will automatically be dropped}. The zero-inflation model uses a logit link.
##' @param dispformula a \emph{one-sided} formula for dispersion containing only fixed effects: the default \code{~1} specifies the standard dispersion given any family. The argument is ignored for families that do not have a dispersion parameter. For an explanation of the dispersion parameter for each family, see \code{\link{sigma}}. The dispersion model uses a log link. In Gaussian mixed models, \code{dispformula=~0} fixes the residual variance to be 0 (actually a small non-zero value), forcing variance into the random effects. The precise value can be controlled via \code{control=glmmTMBControl(zero_dispval=...)}; the default value is \code{sqrt(.Machine$double.eps)}.
##' @param fr model frame
##' @param yobs observed y
##' @param size number of trials in binomial and betabinomial families
##' @param family family object
##' @param start starting values, expressed as a list with possible components \code{beta}, \code{betazi}, \code{betad} (fixed-effect parameters for conditional, zero-inflation, dispersion models); \code{b}, \code{bzi} (conditional modes for conditional and zero-inflation models); \code{theta}, \code{thetazi} (random-effect parameters, on the standard deviation/Cholesky scale, for conditional and z-i models); \code{psi} (extra family parameters, e.g., shape for Tweedie models).
##' @param sparseX see \code{\link{glmmTMB}}
##' @param start_method Options to initialise the starting values for rr parameters; jitter.sd adds variation to the starting values of latent variables when start = "res".
##' @keywords internal
##' @importFrom stats ppois pbinom rnorm
startParams <- function(parameters,
                        formula, ziformula, dispformula,
                        fr,
                        yobs,
                        weights,
                        size = NULL,
                        Xd = NULL,
                        XdS = NULL,
                        family,
                        condReStruc,
                        start = NULL,
                        sparseX = NULL,
                        start_method = list(method = NULL, jitter.sd = 0)) {

  start.met <- start_method$method
  jitter.sd <- ifelse(!is.null(start_method$jitter.sd),
                        start_method$jitter.sd, 0)

    ## rrValues calculates residuals from the fixed model,
    ## fits a reduced rank model to obtain starting values for the latent variables and the factor loadings
    rrValues <- function(yobs, weights, fr, mu,
                       family, formula, ziformula, dispformula, condReStruc,
                       phi = NULL, jitter.sd = 0){
    nobs <- length(yobs)
    resid <- rep(NA, nobs)

    # get either dunn-smyth residuals or
    fam <- family$family
    res.families <- c("poisson", "nbinom2", "binomial", "gaussian")
    if (fam %in% res.families) {
      #### Get the dunn smyth residuals
      if (fam == "poisson") {
        a <- ppois(yobs - 1, mu)
        b <- ppois(yobs, mu)
        u <- runif(n = nobs, min = a, max = b)
        resid <- qnorm(u)
      }
      if (fam == "nbinom2") {
        phi <- phi + 1e-05
        a <- pnbinom(yobs - 1, mu =  mu, size = phi)
        b <- pnbinom(yobs, mu =  mu, size = phi)
        u <- runif(n = nobs, min = a, max = b)
        resid <- qnorm(u)
      }
      if (fam == "nbinom1") {
        phi <- phi + 1e-05
        a <- pnbinom(yobs - 1, mu =  mu, size = mu/phi)
        b <- pnbinom(yobs, mu =  mu, size = mu/phi)
        u <- runif(n = nobs, min = a, max = b)
        resid <- qnorm(u)
      }
      if (fam == "binomial"){
        a <- pbinom(yobs - 1, 1, mu)
        b <- pbinom(yobs, 1, mu)
        u <- runif(n = nobs, min = a, max = b)
        resid <- qnorm(u)
      }
      if (fam == "gaussian"){
        resid <- yobs - mu
      }
    } else {
      resid <- family$dev.resids(y = yobs, mu = mu, wt = weights)
    }
    resid[is.infinite(resid)] <- 0; resid[is.nan(resid)] <- 0
    resid <- as.data.frame(resid)

    get_rank <- function(x) {
      if (x[["blockCode"]] != .valid_covstruct[["rr"]]) return(0)
      p <- x$blockSize
      nt <- x$blockNumTheta
      rank <- (2*p + 1 - sqrt((2*p+1)^2 - 8*nt))/2
      return(rank)
    }

    rank <- vapply(condReStruc,
                   get_rank,
                   FUN.VALUE=numeric(1))
    nlv <- rank[rank > 0]
    namBlk <- names(nlv)

    par.list <- vector("list", length = 3)
    names(par.list) <- c("theta", "b", "fact_load")
    # Use glmmTMB to get initial starting values for factor loadings and latent variables
    fr.res <- cbind(fr, resid)
    ranForm <- no_specials(findbars_x(RHSForm(formula)))
    nrr <- length(namBlk)
    rrTrm <- lapply(1:length(namBlk), function(x) as.character(ranForm[ranForm == namBlk][[x]]))
    x <- sapply(1:nrr, function(x) paste(rrTrm[[x]][2], rrTrm[[x]][1], rrTrm[[x]][3]))
    resForm <- formula(paste("resid ~ 0 "))
    for(i in 1:nrr){
      rrForm <- formula(paste("~ rr(", x[i], ",", nlv[i], ")"))
      resForm <- addForm(resForm, rrForm)
    }
    # residual model; assuming gaussian and fixing sd to 1
    fit.res <- glmmTMB(resForm, data = fr.res, family = gaussian, start = list(betad = c(log(1))), map = list(betad = factor(c(NA))))
    par.list$theta <- fit.res$obj$env$parList(fit.res$fit$par, fit.res$fit$parfull)$theta
    par.list$b <- fit.res$obj$env$parList(fit.res$fit$par, fit.res$fit$parfull)$b
    # Add jitter to latent variables
    par.list$b <- par.list$b + rnorm(length(par.list$b), 0, jitter.sd)

    return(par.list)
  }

  # Fit a fixed model to get the starting values for the fixed parameters, call rrValues to get starting parameters for the rr cov struct (theta and b)
  startVals <-  function(yobs, weights, fr, Xd, XdS, sparseX,
                         family, formula, ziformula, dispformula, condReStruc,
                         parameters, jitter.sd){
    start <- parameters #starting parameters
    fam <- family$family
    ### fit a fixed model
    fixedform <- formula
    RHSForm(fixedform) <- nobars(RHSForm(fixedform))
    # FIX ME: Need to add offset?
    fit.fixed <- glmmTMB(fixedform, data = fr, family = fam,
                         ziformula = ziformula, dispformula = dispformula,
                         weights = weights, sparseX = sparseX)
    fixed.pars <- fit.fixed$obj$env$parList(fit.fixed$fit$par, fit.fixed$fit$parfull)
    nu <- predict(fit.fixed)
    mu <- family$linkinv(nu)

    sparseXd <- ifelse(dim(Xd)[1] == 0 && dim(Xd)[2] == 0, 1, 0)
    if(length(fixed.pars$betad) != 0){
      if(!sparseXd)
        phi <- as.matrix(Xd) %*% exp(fixed.pars$betad)
      else
        phi <- as.vector(XdS %*% exp(fixed.pars$betad))
    }
    # obtain residuals and get starting values for rr
    rrStart <- rrValues(yobs, weights, fr, mu,
                        family, formula, ziformula, dispformula, condReStruc,
                        phi, jitter.sd)

    start.fixed <- fixed.pars
    # Set starting values for fixed parameters from model fit.fixed
    fix.names <- !(names(start) %in% c("b", "theta"))
    for (i in names(start)[fix.names]) {
      if (length(start[[i]]) > 0 & (length(start.fixed[[i]]) == length(start[[i]])))
        start[[i]] <- start.fixed[[i]]
    }

    # Change starting parameters for b and theta for the rr structure
    tp <- trrp <- 1  #theta position for full model, and for model with only rr
    bp <- brrp <- 1  #b position for full model, and for model with only rr
    for (j in seq_along(condReStruc)) {
      nt <- condReStruc[[j]]$blockNumTheta
      nb <- condReStruc[[j]]$blockReps * condReStruc[[j]]$blockSize
      if (condReStruc[[j]]$blockCode == 9) {
        start$b[bp:(bp + nb - 1)] <- rrStart$b[brrp:(brrp + nb - 1)]
        start$theta[tp:(tp + nt - 1)] <- rrStart$theta[trrp:(trrp + nt - 1)]
        brrp <- brrp + nb; trrp <- trrp + nt
      }
      bp <- bp + nb
      tp <- tp + nt
    }

    return(start)
  }

  if(!is.null(start.met)){
    start <- startVals(yobs, weights, fr, Xd, XdS, sparseX,
                       family, formula, ziformula, dispformula, condReStruc,
                       parameters, jitter.sd)
  }

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

  return(parameters)
}

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
##' @param sparseX see \code{\link{glmmTMB}}
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
                       map=NULL,
                       sparseX=NULL,
                       control=glmmTMBControl()) {

  ## handle family specified as naked list
  ## if specified as character or function, should have been converted
  ## to a list with class "family" above ...

  if (is.null(sparseX)) sparseX <- logical(0)
  for (component in c("cond", "zi", "disp")) {
      if (!component %in% names(sparseX)) {
          sparseX[[component]] <- FALSE
      }
  }
  ## make sure the order is right!
  sparseX <- sparseX[c("cond","zi","disp")]

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
        betad_init <- control$zerodisp_val
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

    condList  <- getXReTrms(formula, mf, fr, type="conditional", contrasts=contrasts, sparse=sparseX[["cond"]])
    ziList    <- getXReTrms(ziformula, mf, fr, type="zero-inflation", contrasts=contrasts, sparse=sparseX[["zi"]])
    dispList  <- getXReTrms(dispformula, mf, fr,
                            ranOK=FALSE, type="dispersion",
                            contrasts=contrasts, sparse=sparseX[["disp"]])

    condReStruc <- with(condList, getReStruc(reTrms, ss, aa, reXterms, fr))
    ziReStruc <- with(ziList, getReStruc(reTrms, ss, aa, reXterms, fr))

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


  denseXval <- function(component,lst) if (sparseX[[component]]) matrix(nrow=0,ncol=0) else lst$X
  ## need a 'dgTMatrix' (double, general, Triplet representation)
  sparseXval <- function(component,lst) {
    if (sparseX[[component]]) lst$X else nullSparseMatrix()
  }

  data.tmb <- namedList(
    X = denseXval("cond",condList),
    XS = sparseXval("cond",condList),
    Z = condList$Z,
    Xzi = denseXval("zi",ziList),
    XziS = sparseXval("zi",ziList),
    Zzi = ziList$Z,
    Xd = denseXval("disp",dispList),
    XdS = sparseXval("disp",dispList),

    ## Zdisp=dispList$Z,
    ## use c() on yobs, size to strip attributes such as 'AsIs'
    ##  (which confuse MakeADFun)
    yobs = c(yobs),
    respCol,
    ## strip attributes from offset terms
    offset = c(condList$offset),
    zioffset = c(ziList$offset),
    doffset = c(dispList$offset),
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

  # function to set value for dorr
  rrVal <- function(lst) if(any(lst$ss == "rr")) 1 else 0
  dorr = rrVal(condList)

  getVal <- function(obj, component)
    vapply(obj, function(x) x[[component]], numeric(1))

  ## safer initialization for link functions that might give
  ##  illegal predictions for certain families
  ##  (sqrt() behaves weirdly for beta=0
  ##    [because inverse-link is symmetric around 0?]
  beta_init <-  if (family$link %in% c("identity","inverse","sqrt")) 1 else 0

  ## Extra family specific parameters
  ## FIXME: switch/rewrite to be less ugly?
  numThetaFamily <- if (family$family %in% c("t", "tweedie"))
                    { 1 } else if (family$family == "ordbeta") { 2 } else { 0 }

  rr0 <- function(n) {
       if (is.null(n)) numeric(0) else rep(0, n)
  }

  # theta is 0, except if dorr, theta is 1
  t01 <- function(dorr, condReStruc){
    theta <- rr0(sum(getVal(condReStruc,"blockNumTheta")))
    if(dorr){
      nt <- 1
      blockNumTheta <- getVal(condReStruc,"blockNumTheta")
      blockCode <- getVal(condReStruc, "blockCode")
      for (i in 1:length(blockCode)) {
        if(blockCode[i]==9){
          theta[nt:(nt + blockNumTheta[i] - 1)] <- rep(1, blockNumTheta[i])
        }
        nt <- nt + blockNumTheta[i]
      }
    }
    theta
  }

  parameters <- with(data.tmb,
                     list(
                       beta    = rep(beta_init, max(ncol(X),ncol(XS))),
                       betazi  = rr0(max(ncol(Xzi),ncol(XziS))),
                       b       = rep(beta_init, ncol(Z)),
                       bzi     = rr0(ncol(Zzi)),
                       betad   = rep(betad_init, max(ncol(Xd),ncol(XdS))),
                       theta   = t01(dorr, condReStruc),
                       thetazi = rr0(sum(getVal(ziReStruc,  "blockNumTheta"))),
                       psi  = rr0(numThetaFamily)
                     ))

  if(!is.null(start) || !is.null(control$start_method$method)){
    parameters <- startParams(parameters,
                              formula, ziformula, dispformula,
                              fr,
                              yobs = data.tmb$yobs,
                              weights = data.tmb$weights,
                              size = data.tmb$size,
                              Xd = data.tmb$Xd,
                              XdS = data.tmb$XdS,
                              family,
                              condReStruc,
                              start = start,
                              sparseX = sparseX,
                              start_method = control$start_method)
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
            fr, se, call, verbose, REML, map, sparseX))
}

##' Create X and random effect terms from formula
##' @param formula current formula, containing both fixed & random effects
##' @param mf matched call
##' @param fr full model frame
##' @param ranOK random effects allowed here?
##' @param type label for model type
##' @param contrasts a list of contrasts (see ?glmmTMB)
##' @param sparse (logical) return sparse model matrix?
##' @return a list composed of
##' \item{X}{design matrix for fixed effects}
##' \item{Z}{design matrix for random effects}
##' \item{reTrms}{output from \code{\link{mkReTrms}} from \pkg{lme4}}
##' \item{ss}{splitform of the formula}
##' \item{aa}{additional arguments, used to obtain rank}
##' \item{terms}{terms for the fixed effects}
##' \item{offset}{offset vector, or vector of zeros if offset not specified}
##' \item{reXterms}{terms for the model matrix in each RE term}
##'
##' @importFrom stats model.matrix contrasts
##' @importFrom methods new
##' @importFrom lme4 findbars nobars
getXReTrms <- function(formula, mf, fr, ranOK=TRUE, type="", contrasts, sparse=FALSE) {
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
        terms <- list(fixed=terms(terms_fixed))
        if (!sparse) {
            X <- model.matrix(drop.special(fixedform), fr, contrasts)
        } else {
            X <- Matrix::sparse.model.matrix(drop.special(fixedform), fr, contrasts)
            ## FIXME? ?sparse.model.matrix recommends MatrixModels::model.Matrix(*,sparse=TRUE)
            ##  (but we may not need it, and would add another dependency etc.)
        }
        ## will be 0-column matrix if fixed formula is empty
        offset <- rep(0,nobs)
        if (inForm(fixedform,quote(offset))) {
            ## hate to match offset terms with model frame names
            ##  via deparse, but since that what was presumably done
            ##  internally to get the model frame names in the first place ...
            for (o in extractForm(fixedform,quote(offset))) {
                offset_nm <- deparse1(o)
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

    if (is.null(findbars_x(ranform))) {
        reTrms <- reXterms <- NULL
        Z <- new("dgCMatrix",Dim=c(as.integer(nobs),0L)) ## matrix(0, ncol=0, nrow=nobs)
        aa <- integer(0) #added for rr to get rank
        ss <- integer(0)
    } else {

        ## FIXME: check whether predvars are carried along correctly in terms
        if (!ranOK) stop("no random effects allowed in ", type, " term")
        RHSForm(ranform) <- subbars(RHSForm(reOnly(formula)))

        mf$formula <- ranform
        reTrms <- mkReTrms(no_specials(findbars_x(formula)),
                           fr, reorder.terms=FALSE)

        ss <- splitForm(formula)
        # FIX ME: migrate this (or something like it) down to reTrms,
        ##    allow for more different covstruct types that have additional arguments
        ##  e.g. phylo(.,tree); fixed(.,Sigma)
        # FIX ME: use NA rather than 0 as a placeholder in aa?
        ## FIXME: make sure that eval() happens in the right environment/
        ##    document potential issues
        get_num <- function(v) {
            if (length(v) == 1) return(NA_real_)
            payload <- v[[2]]
            res <- tryCatch(eval(payload, envir = environment(formula)),
                            error = function(e)
                              stop("can't evaluate reduced-rank dimension ",
                                   sQuote(deparse(payload)),
                                   .call = FALSE))
            if (is.na(suppressWarnings(as.numeric(res)))) {
                stop("non-numeric value for reduced-rank dimension",
                     call. = FALSE)
            }
            return(res)
        }
        aa <- ifelse(ss$reTrmClass=="rr",
                     vapply(ss$reTrmAddArgs,
                           get_num,
                           FUN.VALUE=numeric(1)),
                    0)

        ## terms for the model matrix in each RE term
        ## this is imperfect: it should really be done in mkReTrms/mkBlist,
        ## where we are generating these terms anyway on the way
        ## to constructing Z, but that's in lme4 so we can't change it
        ## unless absolutely necessary
        termsfun <- function(x) {
            ## this is a little magic: copying lme4:::mkBlist approach
            ff <- eval(substitute( ~ foo, list(foo = x[[2]]))) ## make formula from LHS
            tt <- try(terms(ff, data=fr), silent=TRUE)         ## construct terms
            if (inherits(tt,"try-error")) {
                stop(
                    sprintf("can't evaluate RE term %s: simplify?",
                            sQuote(deparse(ff)))
                )
            }
            tt
        }

        reXterms <- lapply(ss$reTrmFormulas, termsfun)

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

    namedList(X, Z, reTrms, ss, aa, terms, offset, reXterms)
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
##' @param reXterms terms objects corresponding to each RE term
##' @param fr model frame
##' @param aa additional arguments (i.e. rank)
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
##' @importFrom stats setNames dist .getXlevels
##' @export
getReStruc <- function(reTrms, ss=NULL, aa=NULL, reXterms=NULL, fr=NULL) {

  ## information from ReTrms is contained in cnms, flist elements
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

        if ( any(is.na(aa[ss=="rr"]))) {
          aa0 <- which(is.na(aa) & ss=="rr")
          aa[aa0] <- 2 #set default rank to 2 if it's not specified
        }

        if ( is.null(aa)) {
          aa <- rep(0,length(blksize)) #set rank to 0
        }

        blkrank <- aa
        covCode <- .valid_covstruct[ss]

        parFun <- function(struc, blksize, blkrank) {
            switch(as.character(struc),
                   "0" = blksize, # diag
                   "1" = blksize * (blksize+1) / 2, # us
                   "2" = blksize + 1, # cs
                   "3" = 2,  # ar1
                   "4" = 2,  # ou
                   "5" = 2,  # exp
                   "6" = 2,  # gau
                   "7" = 3,  # mat
                   "8" = 2 * blksize - 1, # toep
                   "9" = blksize * blkrank - (blkrank - 1) * blkrank / 2) #rr
        }
        blockNumTheta <- mapply(parFun, covCode, blksize, blkrank, SIMPLIFY=FALSE)

        ans <- list()
        for (i in seq_along(ss)) {
            tmp <- list(blockReps = nreps[i],
                        blockSize = blksize[i],
                        blockNumTheta = blockNumTheta[[i]],
                        blockCode = covCode[i]
                        )
            if(ss[i] == "ar1") {
                ## FIXME: Keep this warning ?
                if (any(reTrms$cnms[[i]][1] == "(Intercept)") )
                    warning("AR1 not meaningful with intercept")
                if (length(.getXlevels(reXterms[[i]],fr))!=1) {
                    stop("ar1() expects a single, factor variable as the time component")
                }
            } else if(ss[i] == "ou"){
                times <- parseNumLevels(reTrms$cnms[[i]])
                if (ncol(times) != 1)
                    stop("'ou' structure is for 1D coordinates only.")
                if (is.unsorted(times, strictly=TRUE))
                    stop("'ou' is for strictly sorted times only.")
                tmp$times <- drop(times)
            } else if(ss[i] %in% c("exp", "gau", "mat")){
                coords <- parseNumLevels(reTrms$cnms[[i]])
                tmp$dist <- as.matrix( dist(coords) )
            }
            ans[[i]] <- tmp
        }
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
##' @param data data frame (tibbles are OK) containing model variables. Not required, but strongly recommended; if \code{data} is not specified, downstream methods such as prediction with new data (\code{predict(fitted_model, newdata = ...)}) will fail. If it is necessary to call \code{glmmTMB} with model variables taken from the environment rather than from a data frame, specifying \code{data=NULL} will suppress the warning message.
##' @param family a family function, a character string naming a family function, or the result of a call to a family function (variance/link function) information. See \code{\link{family}} for a generic discussion of families or \code{\link{family_glmmTMB}} for details of \code{glmmTMB}-specific families.
##' @param ziformula a \emph{one-sided} (i.e., no response variable) formula for zero-inflation combining fixed and random effects: the default \code{~0} specifies no zero-inflation. Specifying \code{~.} sets the zero-inflation formula identical to the right-hand side of \code{formula} (i.e., the conditional effects formula); terms can also be added or subtracted. \strong{When using \code{~.} as the zero-inflation formula in models where the conditional effects formula contains an offset term, the offset term will automatically be dropped}. The zero-inflation model uses a logit link.
##' @param dispformula a \emph{one-sided} formula for dispersion containing only fixed effects: the default \code{~1} specifies the standard dispersion given any family. The argument is ignored for families that do not have a dispersion parameter. For an explanation of the dispersion parameter for each family, see \code{\link{sigma}}. The dispersion model uses a log link. In Gaussian mixed models, \code{dispformula=~0} fixes the residual variance to be 0 (actually a small non-zero value), forcing variance into the random effects. The precise value can be controlled via \code{control=glmmTMBControl(zero_dispval=...)}; the default value is \code{sqrt(.Machine$double.eps)}.
##' @param weights weights, as in \code{glm}. Not automatically scaled to have sum 1.
##' @param offset offset for conditional model (only).
##' @param contrasts an optional list, e.g., \code{list(fac1="contr.sum")}. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
##' @param na.action a function that specifies how to handle observations
##' containing \code{NA}s.  The default action (\code{na.omit},
##' inherited from the 'factory fresh' value of
##' \code{getOption("na.action")}) strips any observations with any
##' missing values in any variables. Using \code{na.action = na.exclude}
##' will similarly drop observations with missing values while fitting the model,
##' but will fill in \code{NA} values for the predicted and residual
##' values for cases that were excluded during the fitting process
##' because of missingness.
##' @param se whether to return standard errors.
##' @param verbose whether progress indication should be printed to the console.
##' @param doFit whether to fit the full model, or (if FALSE) return the preprocessed data and parameter objects, without fitting the model.
##' @param control control parameters, see \code{\link{glmmTMBControl}}.
##' @param REML whether to use REML estimation rather than maximum likelihood.
##' @param start starting values, expressed as a list with possible components \code{beta}, \code{betazi}, \code{betad} (fixed-effect parameters for conditional, zero-inflation, dispersion models); \code{b}, \code{bzi} (conditional modes for conditional and zero-inflation models); \code{theta}, \code{thetazi} (random-effect parameters, on the standard deviation/Cholesky scale, for conditional and z-i models); \code{psi} (extra family parameters, e.g., shape for Tweedie models).
##' @param map a list specifying which parameter values should be fixed to a constant value rather than estimated. \code{map} should be a named list containing factors corresponding to a subset of the internal parameter names (see \code{start} parameter). Distinct factor values are fitted as separate parameter values, \code{NA} values are held fixed: e.g., \code{map=list(beta=factor(c(1,2,3,NA)))} would fit the first three fixed-effect parameters of the conditional model and fix the fourth parameter to its starting value. In general, users will probably want to use \code{start} to specify non-default starting values for fixed parameters. See \code{\link[TMB]{MakeADFun}} for more details.
##' @param sparseX a named logical vector containing (possibly) elements named "cond", "zi", "disp" to indicate whether fixed-effect model matrices for particular model components should be generated as sparse matrices, e.g. \code{c(cond=TRUE)}. Default is all \code{FALSE}
##' @importFrom stats gaussian binomial poisson nlminb as.formula terms model.weights
##' @importFrom lme4 subbars findbars mkReTrms nobars
##' @importFrom Matrix t
##' @importFrom TMB MakeADFun sdreport
##' @details
##' \itemize{
##' \item Binomial models with more than one trial (i.e., not binary/Bernoulli) can either be specified in the form \code{prob ~ ..., weights = N}, or in the more typical two-column matrix \code{cbind(successes,failures)~...} form.
##' \item Behavior of \code{REML=TRUE} for Gaussian responses matches \code{lme4::lmer}. It may also be useful in some cases with non-Gaussian responses (Millar 2011). Simulations should be done first to verify.
##' \item Because the \code{\link{df.residual}} method for \code{glmmTMB} currently counts the dispersion parameter, users should multiply this value by \code{sqrt(nobs(fit) / (1+df.residual(fit)))} when comparing with \code{lm}.
##' \item Although models can be fitted without specifying a \code{data} argument, its use is strongly recommended; drawing model components from the global environment, or using \code{df$var} notation within model formulae, can lead to confusing (and sometimes hard-to-detect) errors.
##' \item By default, vector-valued random effects are fitted with unstructured (general symmetric positive definite) variance-covariance matrices. Structured variance-covariance matrices can be specified in the form \code{struc(terms|group)}, where \code{struc} is one of
##' \itemize{
##' \item \code{diag} (diagonal, heterogeneous variance)
##' \item \code{ar1} (autoregressive order-1, homogeneous variance)
##' \item \code{cs} (compound symmetric, heterogeneous variance)
##' \item \code{ou} (* Ornstein-Uhlenbeck, homogeneous variance)
##' \item \code{exp} (* exponential autocorrelation)
##' \item \code{gau} (* Gaussian autocorrelation)
##' \item \code{mat} (* Matérn process correlation)
##' \item \code{toep} (* Toeplitz)
##' \item \code{rr} (reduced rank/factor-analytic model)
##' }
##' Structures marked with * are experimental/untested. See \code{vignette("covstruct", package = "glmmTMB")} for more information.
##' \item For backward compatibility, the \code{family} argument can also be specified as a list comprising the name of the distribution and the link function (e.g. \code{list(family="binomial", link="logit")}). However, \strong{this alternative is now deprecated}; it produces a warning and will be removed at some point in the future. Furthermore, certain capabilities such as Pearson residuals or predictions on the data scale will only be possible if components such as \code{variance} and \code{linkfun} are present, see \code{\link{family}}.
##' }
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
##' \donttest{
##' (m1 <- glmmTMB(count ~ mined + (1|site),
##'   zi=~mined,
##'   family=poisson, data=Salamanders))
##' summary(m1)
##' ##' ## Zero-inflated negative binomial model
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
##'
##'
##' ## Using 'map' to fix random-effects SD to 10
##' m1_map <- update(m1, map=list(theta=factor(NA)),
##'                  start=list(theta=log(10)))
##' VarCorr(m1_map)
##' }
glmmTMB <- function(
    formula,
    data = NULL,
    family = gaussian(),
    ziformula = ~0,
    dispformula= ~1,
    weights=NULL,
    offset=NULL,
    contrasts=NULL,
    na.action,
    se=TRUE,
    verbose=FALSE,
    doFit=TRUE,
    control=glmmTMBControl(),
    REML=FALSE,
    start=NULL,
    map=NULL,
    sparseX=NULL
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

    if (inForm(formula, quote(`$`))) {
        warning("use of the ", sQuote("$"), " operator in formulas is not recommended")
    }

    if (missing(data)) {
        warning("use of the ", sQuote("data"), " argument is recommended")
    }

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
    ## need to evaluate offset within environment
    ## how do we figure out where offset exists/whether it has
    ## been prematurely evaluated?
    offsub <- substitute(offset)
    if (is.numeric(offsub)) {
        ## length may cause problems in formula
        data[["..offset"]] <- offset
        offsub <-quote(..offset)
    }
    if (!is.null(eval(offsub,data,
                      enclos=environment(formula)))) {
        formula <- addForm0(formula,makeOp(offsub,
                                           op=quote(offset)))
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
    mf$data <- data ## propagate ..offset modification?

    ## replace . in ziformula with conditional formula, ignoring offset
    if (inForm(ziformula,quote(.))) {
        ziformula <-
            update(RHSForm(drop.special(formula),as.form=TRUE),
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
    ## extract response variable
    ## (name *must* be 'y' to match guts of family()$initialize
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
                   map=map,
                   sparseX=sparseX,
                   control=control)

    ## Allow for adaptive control parameters
    TMBStruc$control <- lapply(control, eval, envir = TMBStruc)

    ## short-circuit
    if (!doFit) return(TMBStruc)

    ## pack all the bits we will need for fitTMB
    res <- fitTMB(TMBStruc)
    return(res)
}

##' Control parameters for glmmTMB optimization
##' @param optCtrl   Passed as argument \code{control} to optimizer. Default value (if default \code{nlminb} optimizer is used): \code{list(iter.max=300, eval.max=400)}
##' @param optArgs   additional arguments to be passed to optimizer function (e.g.: \code{list(method="BFGS")} when \code{optimizer=optim})
##' @param profile   (logical) Experimental option to improve speed and
##'                  robustness when a model has many fixed effects
##' @param collect   (logical) Experimental option to improve speed by
##'                  recognizing duplicated observations.
##' @param parallel  (integer) Set number of OpenMP threads to evaluate
##' the negative log-likelihood in parallel. The default is to evaluate
##' models serially (i.e. single-threaded); users can set a default value
##' for an R session via \code{options(glmmTMB.cores=<value>)}. At present
##' reduced-rank models (i.e., a covariance structure using \code{rr(...)})
##' cannot be fitted in parallel; the number of threads will be automatically
##' set to 1, with a warning if this overrides the user-specified value.
##' @param optimizer Function to use in model fitting. See \code{Details} for required properties of this function.
##' @param eigval_check Check eigenvalues of variance-covariance matrix? (This test may be very slow for models with large numbers of fixed-effect parameters.)
##' @param zerodisp_val value of the dispersion parameter when \code{dispformula=~0} is specified
##' @param start_method (list) Options to initialize the starting values when fitting models with reduced-rank (\code{rr}) covariance structures; \code{jitter.sd} adds variation to the starting values of latent variables when \code{method = "res"}.
##' @param rank_check Check whether all parameters in fixed-effects models are identifiable? This test may be slow for models with large numbers of fixed-effect parameters, therefore default value is 'warning'. Alternatives include 'skip' (no check), 'stop' (throw an error), and 'adjust' (drop redundant columns from the fixed-effect model matrix).
##' @param conv_check Do basic checks of convergence (check for non-positive definite Hessian and non-zero convergence code from optimizer). Default is 'warning'; 'skip' ignores these tests (not recommended for general use!)
##' @details
##' By default, \code{\link{glmmTMB}} uses the nonlinear optimizer
##' \code{\link{nlminb}} for parameter estimation. Users may sometimes
##' need to adjust optimizer settings in order to get models to
##' converge. For instance, the warning \sQuote{iteration limit reached
##' without convergence} may be fixed by increasing the number of
##' iterations using (e.g.)
##'
##' \code{glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3))}.
##'
##' Setting \code{profile=TRUE} allows \code{glmmTMB} to use some special
##' properties of the optimization problem in order to speed up estimation
##' in cases with many fixed effects.
##'
##' Control parameters may depend on the model specification. The value
##' of the controls is evaluated inside an R object that is derived from
##' the output of the \code{\link{mkTMBStruc}} function. For example,
##' to specify that \code{profile} should be enabled if the model has
##' more than 5 fixed-effect parameters, specify
##'
##' \code{profile=quote(length(parameters$beta)>=5)}
##'
##' The \code{optimizer} argument can be any optimization (minimizing) function, provided that:
##' \itemize{
##' \item the first three arguments, in order, are the starting values, objective function, and gradient function;
##' \item the function also takes a \code{control} argument;
##' \item the function returns a list with elements (at least) \code{par}, \code{objective}, \code{convergence} (0 if convergence is successful) and \code{message}
##' (\code{glmmTMB} automatically handles output from \code{optim()}, by renaming the \code{value} component to \code{objective})
##'
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
                           parallel = getOption("glmmTMB.cores", 1L),
                           eigval_check = TRUE,
                           zerodisp_val=log(sqrt(.Machine$double.eps)),
                           start_method = list(method = NULL, jitter.sd = 0),
                           rank_check = c("warning", "adjust", "stop", "skip"),
                           conv_check = c("warning", "skip")) {

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

    rank_check <- match.arg(rank_check)
    conv_check <- match.arg(conv_check)

    ## FIXME: Change defaults - add heuristic to decide if 'profile' is beneficial.
    ##        Something like
    ## profile = (length(parameters$beta) >= 2) &&
    ##           (family$family != "tweedie")
    ## (TMB tweedie derivatives currently slow)
    namedList(optCtrl, profile, collect, parallel, optimizer, optArgs,
              eigval_check, zerodisp_val, start_method, rank_check, conv_check)
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

##' Adjust a model matrix
##' When not rank deficient, do nothing.
##' When rank deficient matrix, drop columns.
##'
##' @param X               model matrix
##' @param tol             non-negative tolerance for testing for "practically zero" singular values (passed to Matrix::rankMatrix())
##' @param why_dropped     logical indicating whether or not to provide information about sets of collinear predictors (not yet implemented)
##'
##' @importFrom Matrix crossprod diag qr qr2rankMatrix
##' @keywords internal
.adjustX <- function(X, tol=NULL, why_dropped=FALSE){
  # perform QR decomposition
  qr_X <- Matrix::qr(X)
  # check if adjustment is necessary
  if(Matrix::qr2rankMatrix(qr_X) < ncol(X)){
    # base qr
    if(is.qr(qr_X)){
      # use column pivoting from base::qr() to identify which columns to keep and which to drop
      to_keep <- qr_X$pivot[1L:qr_X$rank]
      to_drop <- qr_X$pivot[(qr_X$rank+1L):length(qr_X$pivot)]
    } # sparseQR
    else{
      # diagonal elements of R from QR decomposition to identify which columns to keep and which to drop
      R_diag <- abs(Matrix::diag(qr_X@R))
      # borrowing tolerance criterion from that specified in Matrix::qr2rankMatrix
      if(is.null(tol)) tol <- max(dim(X)) * .Machine$double.eps * max(R_diag)
      to_keep <- which(R_diag >= tol)
      to_drop <- which(R_diag < tol)
    }
    # drop columns
    X <- X[,to_keep,drop=FALSE]
    # if(why_dropped){
    # #   TODO: add message describing WHY columns were dropped
    # #   current idea is to follow process outlined at https://stackoverflow.com/a/74103797/20267047 or derivative thereof
    # #   we will eventually use t(X) %*% X to determine collinearity among predictors
    #   tXX <- Matrix::crossprod(X)
    # }
  }
  return(X)
}

##' Check for identifiability of fixed effects matrices X, Xzi, Xd.
##' When rank_check='adjust', drop columns in X and remove associated parameters.
##' @importFrom Matrix rankMatrix
##' @keywords internal
.checkRankX <- function(TMBStruc, rank_check=c('warning','adjust','stop','skip')) {
  rank_check <- match.arg(rank_check)
  Xnames <- c(conditional = "X", conditional = "XS", "zero-inflation" = "Xzi", "zero-inflation" = "XziS", dispersion = "Xd", dispersion = "XdS")
  betanames <- gsub("X", "beta",
                    gsub("S", "", Xnames))
  # use svd-based Matrix::rankMatrix(X) if we wish to abort or warn
  # FIXME: possibly should be an lapply? but I wanted easy access to nm to make error and warnings more informative
  if(rank_check %in% c('stop', 'warning')){
    for (whichX in Xnames) {
      # only attempt rankMatrix if the X matrix contains info
      if(prod(dim(TMBStruc$data.tmb[[whichX]])) == 0) next
        ## if X is rank deficient, stop or throw a warning
        if (Matrix::rankMatrix(TMBStruc$data.tmb[[whichX]]) < ncol(TMBStruc$data.tmb[[whichX]])){
          # determine the model type for a more indicative error or warning message
          model_type <- names(Xnames)[match(whichX, Xnames)]
          action <- get(rank_check, "package:base")
          action("fixed effects in ", model_type," model are rank deficient")
        } ## if rank-deficient
      } ## loop over X components
  } else
    if(rank_check == 'adjust'){
      for(whichX in Xnames){
        # only attempt adjutment if the X matrix contains info
        if(prod(dim(TMBStruc$data.tmb[[whichX]])) == 0) next
        curX <- TMBStruc$data.tmb[[whichX]]
        # use .adjustX to only keep linearly dependent columns of X matrix
        adjX <- .adjustX(curX)
        # if columns were dropped, update variables accordingly
        if(ncol(adjX) != ncol(curX)){
          # inform the user that columns were dropped
          model_type <- names(Xnames)[match(whichX, Xnames)]
          message("dropping columns from rank-deficient ", model_type," model")

          # use colnames of curX and adjX to identify which columns were dropped
          dropped_names   <- setdiff(colnames(curX), colnames(adjX))
          dropped_indices <- match(dropped_names, colnames(curX))

          # retain names of dropped column for use in model output
          attr(adjX, "col.dropped") <- setNames(dropped_indices, dropped_names)

          ## reduce parameters of appropriate component
          beta_name <- betanames[match(whichX, Xnames)]
          kept_indices <- match(colnames(adjX), colnames(curX))
          TMBStruc$parameters[[beta_name]] <- TMBStruc$parameters[[beta_name]][kept_indices]
          TMBStruc$data.tmb[[whichX]] <- adjX
        } ## if matrix was adjusted
      } ## loop over X components
    } ## if rank_check == 'adjust'
  return(TMBStruc)
}

##' Optimize a TMB model and package results
##'
##' This function (called internally by \code{\link{glmmTMB}}) runs
##' the actual model optimization, after all of the appropriate structures
##' have been set up. It can be useful to run \code{\link{glmmTMB}} with
##' \code{doFit=TRUE}, adjust the components as required, and then
##' finish the fitting process with \code{fitTMB} (however, it is the
##' user's responsibility to make sure that any modifications
##' create an internally consistent final fitted object).
##'
##' @param TMBStruc a list contain
##' @examples
##' m0 <- glmmTMB(count ~ mined + (1|site),
##'              family=poisson, data=Salamanders, doFit=FALSE)
##' names(m0)
##' fitTMB(m0)
##' @export
fitTMB <- function(TMBStruc) {

    control <- TMBStruc$control

    has_any_rr <- function(x) {
        any(vapply(x, function(z) z$blockCode == .valid_covstruct[["rr"]],
                   FUN.VALUE = logical(1)))
    }

    if ((has_any_rr(TMBStruc$condReStruc) ||
        has_any_rr(TMBStruc$ziReStruc)) &&
        TMBStruc$control$parallel > 1) {
        warning("rr() not compatible with parallel execution: setting ncores to 1")
        TMBStruc$control$parallel <- 1
    }

    ## Assign OpenMP threads
    n_orig <- openmp(NULL)
    ## Only proceed farther if OpenMP *is* supported ...
    if (n_orig > 0) {
        openmp(n = control$parallel)
        on.exit({
          openmp(n = n_orig)
          })
    }

    if (control $ collect) {
        ## To avoid side-effects (e.g. nobs.glmmTMB), we restore
        ## original data (with duplicates) after fitting.
        data.tmb.old <- TMBStruc$data.tmb
        TMBStruc$data.tmb <- .collectDuplicates(TMBStruc$data.tmb)
    }

    if(control$rank_check %in% c('warning','stop','adjust')){
      TMBStruc <- .checkRankX(TMBStruc, control$rank_check)
    }

    ## avoid repetition; rely on environment for parameters
    optfun <- function() {
        res <- with(obj,
             if( length(par) ) {
                 do.call(control$optimizer,
                         c(list(par, fn, gr,
                                control = control $ optCtrl),
                           control $ optArgs))
             } else {
                 list( par=par, objective=fn(par))
             })
        ## make optim() results look like nlminb() results (which is
        ## what glmmTMB is expecting downstream)
        ## FIXME: what does nloptr output look like?
        ## nlminb components: par, objective, convergence, message
        ## optim components: par, value, counts, convergence, message
        if ("value" %in% names(res)) {
            res$objective <- res$value
            res$value <- NULL
        }
        return(res)
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
          ## pdHess can be FALSE (FIXME: neither of these fallback options is implemented?)
          ##  * Happens for boundary fits (e.g. dispersion close to 0 - see 'spline' example)
          ##    * Option 1: Fall back to old method
          ##    * Option 2: Skip Newton iterations
          for (iter in seq_len(max.newton.steps)) {
            g <- as.numeric( obj$gr(par) )
            if (any(is.na(g)) || max(abs(g)) < newton.tol) break
            par <- par - solve(h, g)
          }
          if (any(is.na(g))) {
            warning("a Newton step failed in profiling")
            par <- oldpar
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
        if (is.na(obj$fn(obj$par))) {
            stop("negative log-likelihood is NaN at starting parameter values")
        }
        if (any(is.na(obj$gr(obj$par)))) {
            stop("some elements of gradient are NaN at starting parameter values")
        }
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
    ## generic complex eigenvalue checker
    e_complex_check <- function(ev, tol=sqrt(.Machine$double.eps)) {
        if (is.complex(ev)) {
            if ((maxim <- max(abs(Im(ev)))) > tol) {
                ev <- Re(ev)
            } else {
                stop(sprintf("detected complex eigenvalues of covariance matrix (max(abs(Im))=%g: try se=FALSE?",
                             maxim))
            }
        }
        return(ev)
    }
    if(!is.null(sdr$pdHess) && control$conv_check != "skip") {
       if(!sdr$pdHess) {
          ## double-check (slower, more accurate hessian)
          env <- environment(obj$fn)
          par <- env$last.par.best
          if (!is.null(rr <- env$random)) {
              par <- par[-rr]
          }
          h <- numDeriv::jacobian(obj$gr, par)
          h <- .5 * (h + t(h))  ## symmetrize
          eigs <- eigen(h)
          ## complex-values check should be unnecessary because we
          ## now symmetrize the hessian, but who knows ... ?
          ev <- e_complex_check(eigs$values)
          if (min(ev)>.Machine$double.eps) {
              ## apparently fit is OK after all ...
              sdr$pdHess <- TRUE
              Vtheta <- try(solve(h), silent=TRUE)
              if (!inherits(Vtheta,"try-error")) sdr$cov.fixed[] <- Vtheta
          } else {
              warning(paste0("Model convergence problem; ",
                             "non-positive-definite Hessian matrix. ",
                             "See vignette('troubleshooting')"))
          }
      } else if (control$eigval_check && length(sdr$cov.fixed)>0) {
          eigval <- try(1/eigen(sdr$cov.fixed)$values, silent=TRUE)
          if( is(eigval, "try-error") || ( min(e_complex_check(eigval)) < .Machine$double.eps*10 ) ) {
              warning(paste0("Model convergence problem; ",
                             "extreme or very small eigenvalues detected. ",
                             "See vignette('troubleshooting')"))
          } ## bad eigval
      } ## do eigval check
    } ## pdHess exists

    if ( !is.null(fit$convergence) && fit$convergence != 0 && control$conv_check != "skip") {
        warning("Model convergence problem; ",
                fit$message, ". ",
                "See vignette('troubleshooting')")
    }
    
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
                                map,
                                sparseX,
                                parallel = control$parallel))
    ## FIXME: are we including obj and frame or not?
    ##  may want model= argument as in lm() to exclude big stuff from the fit
    ## If we don't include obj we need to get the basic info out
    ##    and provide a way to regenerate it as necessary
    ## If we don't include frame, then we may have difficulty
    ##    with predict() in its current form

    ## FIXME (rank_check): ret needs to know about dropped columns to use
    ##  them in the summary, etc. At this point in the code, this information
    ##  only remains in condList, ziList, and dispList.

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

    mkCoeftab <- function(coefs,vcovs) {
        p <- length(coefs)
        coefs <- cbind("Estimate" = coefs,
                       "Std. Error" = sqrt(diag(vcovs)))
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

    # FIXME (rank_check): to include dropped predictors in the output, fixef
    #  needs to be able to find out about them

    ff <- fixef(object)
    vv <- vcov(object,include_mapped=TRUE)
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


