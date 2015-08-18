## Copy-paste from cpp-source:
## FIXME: to be read from C++ internals
.valid_family <- c(poisson_family           = 0,
                   binomial_family          = 1,
                   negative_binomial_family = 2,
                   Gamma_family             = 3,
                   beta_family              = 4,
                   gaussian_family          = 5,
                   truncated_poisson_family = 6,
                   trunc_NB_family          = 7,
                   logistic_family          = 8,
                   betabinomial_family      = 9 )
.valid_link <- c(log_link                 = 0,
                 logit_link               = 1,
                 probit_link              = 2,
                 inverse_link             = 3,
                 cloglog_link             = 4,
                 identity_link            = 5 )
names(.valid_family) <- sub("_family","",names(.valid_family))
names(.valid_link) <- sub("_link","",names(.valid_link))

## TODO: Steamline construction of design matrices
## - Can we use lme4 utilities ?

##' main TMB function
##' @param formula combined fixed and random effects formula, following lme4 syntac
##' @param data data frame
##' @param family \code{\link{family}}
##' @param ziformula combined fixed and random effects formula for zero-inflation: the default \code{~0} specifies no zero-inflation
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
    weights,
    offset
    )
{

    ## edited copy-paste from glFormula
    ## glFormula <- function(formula, data=NULL, family = gaussian,
    ##    subset, weights, na.action, offset,
    ##     contrasts = NULL, mustart, etastart,
    ##                      control = glmerControl(), ...) {

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
    fr.form <- subbars(formula) # substitute "|" by "+"
    environment(fr.form) <- environment(formula)
    ## model.frame.default looks for these objects in the environment
    ## of the *formula* (see 'extras', which is anything passed in ...),
    ## so they have to be put there ...
    for (i in c("weights", "offset")) {
        if (!eval(bquote(missing(x=.(i)))))
            assign(i, get(i, parent.frame()), environment(fr.form))
    }
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    ## FIXME: throw an error *or* convert character to factor
    ## convert character vectors to factor (defensive)
    ## fr <- factorize(fr.form, fr, char.only = TRUE)
    ## store full, original formula & offset
    attr(fr,"formula") <- formula
    attr(fr,"offset") <- mf$offset
    n <- nrow(fr)
    ## random effects and terms modules
    reTrms <- mkReTrms(findbars(RHSForm(formula)), fr)
    ## sanity checks (skipped!)
    ## wmsgNlev <- checkNlevels(reTrms$ flist, n = n, control, allow.n = TRUE)
    ## wmsgZdims <- checkZdims(reTrms$Ztlist, n = n, control, allow.n = TRUE)
    ## wmsgZrank <- checkZrank(reTrms$ Zt, n = n, control, nonSmall = 1e6, allow.n = TRUE)

    ## fixed-effects model matrix X - remove random effect parts from formula:
    fixedform <- formula
    RHSForm(fixedform) <- nobars(RHSForm(fixedform))
    mf$formula <- fixedform
    ## re-evaluate model frame to extract predvars component
    fixedfr <- eval(mf, parent.frame())
    attr(attr(fr,"terms"),"predvars.fixed") <-
        attr(attr(fixedfr,"terms"),"predvars")

    ## ran-effects model frame (for predvars)
    ## important to COPY formula (and its environment)?
    ranform <- formula
    RHSForm(ranform) <- subbars(RHSForm(reOnly(formula)))
    mf$formula <- ranform
    ranfr <- eval(mf, parent.frame())
    attr(attr(fr,"terms"), "predvars.random") <-
        attr(terms(ranfr), "predvars")

    ## FIXME: make model matrix sparse?? i.e. Matrix:::sparse.model.matrix(...)
    X <- model.matrix(fixedform, fr, contrasts)

    ## if(is.null(rankX.chk <- control[["check.rankX"]]))
    ## rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    ## X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
    ## if(is.null(scaleX.chk <- control[["check.scaleX"]]))
    ##     scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    ## X <- checkScaleX(X, kind=scaleX.chk)

    ## list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula,
    ##  wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))

    ## ... done with

    ## extract family and link information from family object
    link <- family$link
    family <- family$family   ## overwrites family: original info lost

    ## extract response variable
    yobs <- fr[,attr(fr,"response")]

    ## structure (0 =  diag, 1 = us, 2 = comp symm ...)
    ## no. of REs (length of b vector)
    ## no. of parameters  (length of theta vector)
    ##

    ## FIXME: go back and process zero-inflation formula, if any
    nobs <- nrow(fr)
    Xzi <- matrix(ncol=0,nrow=nobs)
    Zzi <- matrix(ncol=0,nrow=nobs)

    ## Get info on sizes of RE components
    theta <- numeric(length(reTrms$theta))

    ## FAKE EXAMPLE
    ## dd <- expand.grid(year=1:10,site=1:20)
    ## dd$size <- dd$y <- rnorm(nrow(dd))
    ## L <- lFormula(y~(1|site)+(size|year),data=dd)
    ## rr <- L$reTrms

    ## hack names of Ztlist to extract grouping variable of each RE term
    grpVar <- gsub("^[^|]*\\| ","",names(reTrms$Ztlist)) ## remove 1st term+|
    getLevs <- function(i) with(reTrms,
                                length(levels(flist[[grpVar[i]]])))
    nlevs <- sapply(seq_along(reTrms$cnms),getLevs)
    blksize <- sapply(reTrms$Ztlist,nrow)/nlevs
    ## figure out number of parameters from block size + structure type

    ## for now *all* RE are unstructured
    struc <- rep(1,length(nlevs))

    parFun <- function(struc,blksize) {
        switch(as.character(struc),
               "0"=blksize, ## diag
               "1"=blksize*(blksize+1)/2, ## us
               "2"=blksize+1)  ## cs
    }
    npars <- mapply(parFun,struc,blksize)

    Z <- t(reTrms$Zt)

    data.tmb <- namedList(
        X,
        Z,
        Xzi,
        Zzi,
        yobs,
        npars,   ## number of variance-covariance params per term
        nlevs,   ## " " levels " "
        struc,   ## structure code
        blksize, ## block size
        family = .valid_family[family],
        link = .valid_link[link]
        )
    parameters <- list(
        beta     = rep(0, ncol(X)) ,
        b        = rep(0, ncol(Z)) ,
        betazi   = rep(0, ncol(Xzi)),
        bzi      = rep(0, ncol(Zzi)),
        theta    = rep(0, length(reTrms$theta))
        )

    ## FIXME: something about dispersion/model matrix?

    obj <- MakeADFun(data.tmb,
                     parameters,
                     random = c("b","bzi"),
                     profile = NULL, ## TODO: Optionally "beta"
                     silent = FALSE, ## TODO: set to TRUE
                     DLL="glmmTMB")
    obj ## For now give the object without optimizing
}


