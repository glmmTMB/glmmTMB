##' @title Create X and random effect terms from formula
##' @param formula current formula, containing both fixed & random effects
##' @param mf matched call
##' @param fr full model frame
##' @param ranOK random effects allowed here?
##' @param type label for model type
##' @return a list composed of
##' \item{X}{design matrix for fixed effects}
##' \item{Z}{design matrix for random effects}
##' \item{fixedfr}{model frame for fixed effects. Reading this in prevents the predictor function from recalculating bases for splines and orthogonal polynomials. Might have to be removed later}
##' \item{ranfr}{as fixedfr but for random effects}
##' \item{reTrms}{output from mkReTerms from LME4}
getXReTrms <- function(formula, mf, fr, ranOK=TRUE, type="") {
    ## fixed-effects model matrix X -
    ## remove random effect parts from formula:
    fixedform <- formula
    RHSForm(fixedform) <- nobars(RHSForm(fixedform))

    nobs <- nrow(fr)
    ## check for empty fixed form

    if (identical(RHSForm(fixedform),~0) ||
        identical(RHSForm(fixedform),~-1)) {
        X <- NULL
        fixedfr <- NULL
    } else {
        mf$formula <- fixedform
        ## re-evaluate model frame to extract predvars component
        ## in *grandparent* environment
        fixedfr <- eval(mf, parent.frame(2))
        attr(attr(fr,"terms"), "predvars.fixed") <-
            attr(attr(fixedfr,"terms"), "predvars")

        ## FIXME: make model matrix sparse?? i.e. Matrix:::sparse.model.matrix(...)
        X <- model.matrix(fixedform, fr, contrasts)
        ## will be 0-column matrix if fixed formula is empty
    }
    ## ran-effects model frame (for predvars)
    ## important to COPY formula (and its environment)?
    ranform <- formula
    if (is.null(findbars(ranform))) {
        ranfr <- NULL
        reTrms <- NULL
        Z <- matrix(0, ncol=0, nrow=nobs)
        ss <- integer(0)
    } else {
        if (!ranOK) stop("no random effects allowed in ", type, " term")
        RHSForm(ranform) <- subbars(RHSForm(reOnly(formula)))

        mf$formula <- ranform
        ranfr <- eval(mf, parent.frame(2))
        attr(attr(fr,"terms"), "predvars.random") <-
            attr(terms(ranfr), "predvars")
        reTrms <- mkReTrms(findbars(RHSForm(formula)), fr)

        ss <- splitForm(formula)
        ss <- unlist(ss$reTrmClasses)

        Z <- as.matrix(t(reTrms$Zt))
    }

    ## if(is.null(rankX.chk <- control[["check.rankX"]]))
    ## rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    ## X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
    ## if(is.null(scaleX.chk <- control[["check.scaleX"]]))
    ##     scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    ## X <- checkScaleX(X, kind=scaleX.chk)

    ## list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula,
    ##      wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))

    ## FIXME: come back and figure out how we need to store fixedfr and ranfr

    return(namedList(X, Z, fixedfr, ranfr, reTrms, ss))
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Calculate random effect structure
##' calculates number of random effects, number of parameters,
##' blocksize and number of blocks.
##' @param reTrms random-effects terms list
##' @return a list
##' \item{blockNumTheta}{number of variance covariance parameters per term}
##' \item{blockSize}{size (dimension) of one block}
##' \item{blockReps}{number of times the blocks are repeated (levels)}
##' \item{covCode}{structure code}
##' data(sleepstudy,package="lme4")
##' rt <- lme4::lFormula(Reaction~Days+(1|Subject)+(0+Days|Subject),
##'                     sleepstudy)$reTrms
##' rt2 <- lme4::lFormula(Reaction~Days+(Days|Subject),
##'                     sleepstudy)$reTrms
##' getReStruc(rt)                    
getReStruc <- function(reTrms,ss) {

    if (is.null(reTrms)) {
        ReStrucList <- list()
    } else {
        ## Get info on sizes of RE components

        ## hack names of Ztlist to extract grouping variable of each RE term
        ## remove 1st term+|
        grpVar <- gsub("^[^|]*\\| ", "", names(reTrms$Ztlist))
        getLevs <- function(i) with(reTrms, length(levels(flist[[grpVar[i]]])))
        nreps <- sapply(seq_along(reTrms$cnms), getLevs)
        blksize <- sapply(reTrms$Ztlist,nrow) / nreps
        ## figure out number of parameters from block size + structure type

        covCode <- .valid_covstruct[ss]

        parFun <- function(struc, blksize) {
            switch(as.character(struc),
                   "0" = blksize, # diag
                   "1" = blksize * (blksize+1) / 2, # us
                   "2" = blksize + 1) # cs
        }
        blockNumTheta <- mapply(parFun, covCode, blksize, SIMPLIFY=FALSE)

        ReStrucList <- mapply(list,
                              blockReps=nreps,
                              blockSize=blksize,
                              blockNumTheta=blockNumTheta,
                              blockCode=covCode,SIMPLIFY=FALSE)
        names(ReStrucList) <- names(reTrms$Ztlist)
    }

    return(ReStrucList)
}

usesDispersion <- function(x) {
  !x %in% c("binomial","poisson","truncated_poisson")
}

##' @title main TMB function
##' @param formula combined fixed and random effects formula, following lme4 syntac
##' @param data data frame
##' @param family \code{\link{family}}
##' @param ziformula combined fixed and random effects formula for zero-inflation: the default \code{~0} specifies no zero-inflation
##' @param dispformula combined fixed and random effects formula for dispersion: the default \code{~0} specifies no zero-inflation
##' @param weights
##' @param offset
##' @param se whether to return standard errors
##' @param debug whether to return the preprocessed data and parameter objects, without fitting the model
##' @importFrom lme4 subbars findbars mkReTrms nobars
##' @importFrom Matrix t
##' @importFrom TMB MakeADFun sdreport
##' @importMethodsFrom TMB print.sdreport
##' @export
##' @examples
##' data(sleepstudy, package="lme4")
##' glmmTMB(Reaction ~ Days + (1|Subject), sleepstudy, debug=TRUE)
##' glmmTMB(Reaction ~ Days + us(1|Subject), sleepstudy, debug=TRUE)
##' glmmTMB(Reaction ~ Days + diag(1|Subject), sleepstudy, debug=TRUE)
glmmTMB <- function (
    formula,
    data = NULL,
    family = gaussian(),
    ziformula = ~0,
    dispformula= NULL, 
    weights=NULL,
    offset=NULL,
    se=FALSE,
    debug=FALSE
    )
{

    ## edited copy-paste from glFormula
    ## glFormula <- function(formula, data=NULL, family = gaussian,
    ##                       subset, weights, na.action, offset,
    ##                       contrasts = NULL, mustart, etastart,
    ##                       control = glmerControl(), ...) {

    ## FIXME: check for offsets in ziformula/dispformula, throw an error

    call <- mf <- mc <- match.call()
    ## extract family, call lmer for gaussian

    ## FIXME: jump through the usual hoops to allow
    ## character, function, family-object 
    if (grepl("^quasi", family$family))
        stop('"quasi" families cannot be used in glmmtmb')

    ## extract family and link information from family object
    link <- family$link
    family <- family$family # overwrites family: original info lost
    
    if (is.null(dispformula)) {
      dispformula <- if (usesDispersion(family)) ~1 else ~0
    }
    
    ## lme4 function for warning about unused arguments in ...
    ## ignoreArgs <- c("start","verbose","devFunOnly",
    ##   "optimizer", "control", "nAGQ")
    ## l... <- list(...)
    ## l... <- l...[!names(l...) %in% ignoreArgs]
    ## do.call(checkArgs, c(list("glmer"), l...))

    # substitute evaluated version
    ## FIXME: denv leftover from lme4, not defined yet
    
    mc$formula <- formula <- as.formula(formula, env = denv)
    
    
    ## now work on evaluating model frame
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")

   
    ## want the model frame to contain the union of all variables
    ## used in any of the terms
    ## combine all formulas

    formList <- list(formula[[3]], ziformula, dispformula)
    formList <- lapply(formList,
                   function(x) noSpecials(subbars(x),delete=FALSE))
                       ## substitute "|" by "+"; drop special
    formList <- gsub("~", "\\+", lapply(formList,safeDeparse)) # character
    combForm <- reformulate(Reduce(paste,formList),
                            response=deparse(formula[[2]]))
    environment(combForm) <- environment(formula)
    ## model.frame.default looks for these objects in the environment
    ## of the *formula* (see 'extras', which is anything passed in ...),
    ## so they have to be put there ...
    for (i in c("weights", "offset")) {
        if (!eval(bquote(missing(x=.(i)))))
            assign(i, get(i, parent.frame()), environment(combForm))
    }

    mf$formula <- combForm
    fr <- eval(mf, parent.frame())
    ## FIXME: throw an error *or* convert character to factor
    ## convert character vectors to factor (defensive)
    ## fr <- factorize(fr.form, fr, char.only = TRUE)
    ## store full, original formula & offset
    attr(fr,"formula") <- combForm
    attr(fr,"offset") <- mf$offset
    nobs <- nrow(fr)

    condList <- getXReTrms(formula, mf, fr)
    ziList    <- getXReTrms(ziformula, mf, fr)
    dispList  <- getXReTrms(dispformula, mf, fr, ranOK=FALSE, "dispersion")

    ## sanity checks (skipped!)
    ## wmsgNlev <- checkNlevels(reTrms$ flist, n=n, control, allow.n=TRUE)
    ## wmsgZdims <- checkZdims(reTrms$Ztlist, n=n, control, allow.n=TRUE)
    ## wmsgZrank <- checkZrank(reTrms$Zt, n=n, control, nonSmall=1e6, allow.n=TRUE)
    
    ## store info on location of response variable
    respCol <- attr(terms(fr),"response")
    names(respCol) <- names(fr)[respCol]
    
    ## extract response variable
    yobs <- fr[,respCol]

    condReStruc <- with(condList,getReStruc(reTrms,ss))
    ziReStruc <- with(ziList,getReStruc(reTrms,ss))

    ## FIXME: deal with offset in formula
    if (grepl("offset",safeDeparse(formula)))
        stop("offsets within formulas not implemented")
    if (is.null(offset)) offset <- rep(0,nobs))
    
    if (is.null(weights <- fr[["(weights)"]]))
        weights <- rep(1,nobs)
    
    data.tmb <- namedList(
        X = condList$X,
        Z = condList$Z,
        Xzi = ziList$X,
        Zzi = ziList$Z,
        Xd = dispList$X,
        ## Zdisp=dispList$Z,
        yobs,
        offset,
        weights,
        ## information about random effects structure
        terms = condReStruc,
        termszi = ziReStruc,
        family = .valid_family[family],
        link = .valid_link[link]
        )
    getVal <- function(object,component) {
        vapply(object,function(x) x[[component]],numeric(1))
    }
    parameters <- with(data.tmb,
      list(
          beta    = rep(0, ncol(X)),
          b       = rep(0, ncol(Z)),
          betazi  = rep(0, ncol(Xzi)),
          bzi     = rep(0, ncol(Zzi)),
          theta   = rep(0, sum(getVal(condReStruc,"blockNumTheta"))),
          thetazi = rep(0, sum(getVal(ziReStruc,"blockNumTheta"))),
          betad   = rep(0, ncol(Xd))
          ))

    ## short-circuit
    if(debug) return(namedList(data.tmb,parameters))

    randomArg <- NULL
    if (ncol(data.tmb$Z) > 0) randomArg <- c(randomArg,"b")
    if (ncol(data.tmb$Zzi) > 0) randomArg <- c(randomArg,"bzi")

    obj <- MakeADFun(data.tmb,
                     parameters,
                     random = randomArg,
                     profile = NULL, # TODO: Optionally "beta"
                     silent = FALSE, # TODO: set to TRUE
                     DLL = "glmmTMB")
    optTime <- system.time(fit <- with(obj, nlminb(start=par,objective=fn,
                                                   gradient=gr)))
    sdr <- if (se) sdreport(obj) else NULL

    modelInfo <- namedList(nobs,respCol,
         restruc=namedList(condReStruc,ziReStruc),
         allForm=namedList(combForm,formula,
                           ziformula,dispformula))
    ## FIXME: are we including obj and frame or not?  
    ##  may want model= argument as in lm() to exclude big stuff from the fit
    ## If we don't include obj we need to get the basic info out
    ##    and provide a way to regenerate it as necessary
    ## If we don't include frame, then we may have difficulty
    ##    with predict() in its current form
    output <- namedList(obj, fit, sdr, call, frame=fr)
    class(output) <- "glmmTMB"

    return(output)
}
