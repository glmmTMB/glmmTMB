## Copy-paste from cpp-source:
## FIXME: to be read from C++ internals
.valid_family <- c(gaussian_family           = 0,
                   binomial_family           = 100,
                   betabinomial_family       = 101,
                   beta_family               = 200,
                   gamma_family              = 300,
                   poisson_family            = 400,
                   truncated_poisson_family  = 401,
                   nbinom1_family            = 500,
                   nbinom2_family            = 501,
                   truncated_nbinom1_family  = 502,
                   truncated_nbinom2_family  = 503,
                   t_family                  = 600,
                   tweedie_family            = 700 )
.valid_link <- c(log_link                 = 0,
                 logit_link               = 1,
                 probit_link              = 2,
                 inverse_link             = 3,
                 cloglog_link             = 4,
                 identity_link            = 5 )
names(.valid_family) <- sub("_family","",names(.valid_family))
names(.valid_link) <- sub("_link","",names(.valid_link))

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param formula current formula, containing both fixed & random effects
##' @param mf matched call
##' @param fr full model frame
##' @param ranOK random effects allowed here?
##' @param type label for model type
##' @return 
getXReTrms <- function(formula,mf,fr,ranOK=TRUE,type="") {
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
        attr(attr(fr,"terms"),"predvars.fixed") <-
            attr(attr(fixedfr,"terms"),"predvars")

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
        Z <- matrix(0,ncol=0,nrow=nobs)
    } else {
        if (!ranOK) stop("no random effects allowed in ",type," term")
        RHSForm(ranform) <- subbars(RHSForm(reOnly(formula)))

        mf$formula <- ranform
        ranfr <- eval(mf, parent.frame(2)) 
        attr(attr(fr,"terms"), "predvars.random") <-
            attr(terms(ranfr), "predvars")
        reTrms <- mkReTrms(findbars(RHSForm(formula)), fr)

        Z <- as.matrix(t(reTrms$Zt))
    }

    ## if(is.null(rankX.chk <- control[["check.rankX"]]))
    ## rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    ## X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
    ## if(is.null(scaleX.chk <- control[["check.scaleX"]]))
    ##     scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    ## X <- checkScaleX(X, kind=scaleX.chk)

    ## list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula,
    ##  wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))

    ## FIXME: come back and figure out how we need to store fixedfr and ranfr

    return(namedList(X,Z,fixedfr,ranfr,reTrms))
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param reTrms  random-effects terms list
##' @return 
getReStruc <- function(reTrms) {

    if (is.null(reTrms)) {
        blksize <- nreps <- covCode <- blockNumTheta <- integer(0)
    } else {
        
        ## Get info on sizes of RE components
        theta <- numeric(length(reTrms$theta))

        ## hack names of Ztlist to extract grouping variable of each RE term
        grpVar <- gsub("^[^|]*\\| ","",names(reTrms$Ztlist)) ## remove 1st term+|
        getLevs <- function(i) with(reTrms,
                                    length(levels(flist[[grpVar[i]]])))
        nreps <- sapply(seq_along(reTrms$cnms),getLevs)
        blksize <- sapply(reTrms$Ztlist,nrow)/nreps
        ## figure out number of parameters from block size + structure type

        ## for now *all* RE are diagonal
        covCode <- rep(0,length(nreps))

        parFun <- function(struc,blksize) {
            switch(as.character(struc),
                   "0"=blksize, ## diag
                   "1"=blksize*(blksize+1)/2, ## us
                   "2"=blksize+1)  ## cs
        }
        blockNumTheta <- mapply(parFun,covCode,blksize)
    }
    
    return(namedList(blockNumTheta,
                     blockSize=blksize,
                     blockReps=nreps,
                     covCode))
}

##' main TMB function
##' @param formula combined fixed and random effects formula, following lme4 syntac
##' @param data data frame
##' @param family \code{\link{family}}
##' @param ziformula combined fixed and random effects formula for zero-inflation: the default \code{~0} specifies no zero-inflation
##' @param dispformula combined fixed and random effects formula for dispersion: the default \code{~0} specifies no zero-inflation
##' @param weights 
##' @param offset 
##' @importFrom lme4 subbars findbars mkReTrms nobars
##' @importFrom Matrix t
##' @importFrom TMB MakeADFun
##' @export
##' @examples
##' data(sleepstudy,package="lme4")
##' glmmTMB(Reaction~Days+(1|Subject),sleepstudy)
glmmTMB <- function (
    formula,
    data = NULL,
    family = gaussian(),
    ziformula = ~0,
    dispformula= ~1, #FIXME: set appropriate family-specific defaults
    weights,
    offset,
    debug=FALSE
    )
{

    ## edited copy-paste from glFormula
    ## glFormula <- function(formula, data=NULL, family = gaussian,
    ##    subset, weights, na.action, offset,
    ##     contrasts = NULL, mustart, etastart,
    ##                      control = glmerControl(), ...) {

    ## FIXME: check for offsets in ziformula/dispformula, throw an error
    ##
    mf <- mc <- match.call()
    ## extract family, call lmer for gaussian

    if (grepl("^quasi",family$family))
        stop('"quasi" families cannot be used in glmmtmb')

    ## ignoreArgs <- c("start","verbose","devFunOnly",
    ##   "optimizer", "control", "nAGQ")
    ## l... <- list(...)
    ## l... <- l...[!names(l...) %in% ignoreArgs]
    ## do.call(checkArgs, c(list("glmer"), l...))

    mc$formula <- formula <- as.formula(formula, env = denv)    ## substitute evaluated version

    ## now work on evaluating model frame
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")


    ## want the model frame to contain the union of all variables
    ## used in any of the terms
    ## combine all formulas

    formList <- list(formula[[3]],ziformula,dispformula)
    formList <- lapply(formList,subbars) # substitute "|" by "+"
    formList <- gsub("~","\\+",lapply(formList,safeDeparse))  ## character
    combForm <- reformulate(Reduce(paste,formList), 
                            response= deparse(formula[[2]]))
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
    n <- nrow(fr)

    fixedList <- getXReTrms(formula,mf,fr)
    ziList    <- getXReTrms(ziformula,mf,fr)
    dispList  <- getXReTrms(dispformula,mf,fr,ranOK=FALSE,"dispersion")

    ## sanity checks (skipped!)
    ## wmsgNlev <- checkNlevels(reTrms$ flist, n = n, control, allow.n = TRUE)
    ## wmsgZdims <- checkZdims(reTrms$Ztlist, n = n, control, allow.n = TRUE)
    ## wmsgZrank <- checkZrank(reTrms$ Zt, n = n, control, nonSmall = 1e6, allow.n = TRUE)

    ## extract family and link information from family object
    link <- family$link
    family <- family$family   ## overwrites family: original info lost

    ## extract response variable
    yobs <- fr[,attr(terms(fr),"response")]

    fixedReStruc <- getReStruc(fixedList$reTrms)
    ziReStruc <- getReStruc(ziList$reTrms)

    data.tmb <- namedList(
        X=fixedList$X,
        Z=fixedList$Z,
        Xzi=ziList$X,
        Zzi=ziList$Z,
        Xd=dispList$X,
        ## Zdisp=dispList$Z,
        yobs,
        ## offset,

        ## information about random effects structure
        blockReps=fixedReStruc$blockReps,     ## nreps,   ## " " levels " "
        blockSize=fixedReStruc$blockSize,     ## blksize, ## block size
        blockNumTheta=fixedReStruc$blockNumTheta, ##  number of variance-covariance params per term
        blockCode=fixedReStruc$covCode,       ## struc,   ## structure code

        blockRepszi=ziReStruc$blockReps,     ## nreps,   ## " " levels " "
        blockSizezi=ziReStruc$blockSize,     ## blksize, ## block size
        ## FIXME: change blockNumTheta to numTheta???
        blockNumThetazi=ziReStruc$blockNumTheta, ##  number of variance-covariance params per term
        blockCodezi=ziReStruc$covCode,       ## struc,   ## structure code



        family = .valid_family[family],
        link = .valid_link[link]
        )
    parameters <- with(data.tmb,
      list(
          beta     = rep(0, ncol(X)) ,
          b        = rep(0, ncol(Z)) ,
          betazi   = rep(0, ncol(Xzi)),
          bzi      = rep(0, ncol(Zzi)),
          theta    = rep(0, sum(blockNumTheta)),
          thetazi  = rep(0, sum(blockNumThetazi)),
          betad    = rep(0, ncol(Xd))
          )
                       )

    ## short-circuit
    if(debug) return(namedList(data.tmb,parameters))

    obj <- MakeADFun(data.tmb,
                     parameters,
                     random = c("b","bzi"),
                     profile = NULL, ## TODO: Optionally "beta"
                     silent = FALSE, ## TODO: set to TRUE
                     DLL="glmmTMB")

    obj ## For now give the object without optimizing

    optTime <- system.time(fit <- with(obj,nlminb(start=par,objective=fn,
                                                  gradient=gr)))
    #sdr <- sdreport(obj)
	
#    return(list(optTime, fit, sdr))
    return(namedList(optTime, fit))
    ## now structure the output object

}


