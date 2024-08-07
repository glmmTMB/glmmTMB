##' Hat values (and other influence measures) for glmmTMB objects
##' 
##' @param model a \code{glmmTMB} model
##' @param diag Get diagonal only?
##' @param \dots unused: for method compatibility only
##' @examples
##' fm0 <- glmmTMB(weight ~ Time, data=ChickWeight)
##' fm1 <- lm(weight ~ Time, data = ChickWeight)
##' all.equal(hatvalues(fm0), unname(hatvalues(fm1)))
##' fm2 <- glmmTMB(weight ~ Time + diag(Time | Chick), data=ChickWeight)
##' fm3 <- lme4::lmer(weight ~ Time + (Time || Chick), data = ChickWeight)
##' all.equal(hatvalues(fm2), unname(hatvalues(fm3)), tolerance = 1e-2)
##' @export
hatvalues.glmmTMB <- function(model, diag=TRUE, ...) {

    ## FIXME: don't depend on RTMB
    ##  (GetTape, MakeTape, TMB::config() below)
    if (requireNamespace("RTMB")) {
        dyn.load(system.file("libs", TMB::dynlib("RTMB"), package = "RTMB"))
    } else {
        stop("sorry, hatvalues currently require that the RTMB package be installed")
    }
    
    check_dots(...)
    has.random <- any(model$obj$env$lrandom())
    obj <- model$obj
    ## We mess with these... (cleanup on exit!)
    restore.on.exit <- c("ADreport",
                         "parameters",
                         "data")
    oldvars <- sapply(restore.on.exit, get, envir=obj$env, simplify=FALSE)
    restore.oldvars <- function(){
        for(var in names(oldvars)) assign(var, oldvars[[var]], envir=obj$env)
    }
    on.exit({restore.oldvars(); obj$retape()})
    ## #################################################################
    ## Get derivatives of prediction
    ##
    ##    mu_hat( b_hat( theta_hat(yobs), yobs) , theta_hat(yobs) )
    ##
    ## wrt. yobs.
    ## Note the three partial derivative 'paths' to consider.
    ## We can get this derivative by
    ##  1. yobs -> theta_hat(yobs)
    ##  2. theta -> mu_hat( b_hat( theta, yobs) , theta ) [fixed yobs]
    ##  3. yobs -> mu_hat( b_hat( theta, yobs) , theta ) [fixed theta]
    ## #################################################################
    ##parhat <- obj$env$last.par.best
    parhat <- model$fit$parfull
    pl <- obj$env$parList(par=parhat)
    yobs <- obj$env$data$yobs
    obj$env$parameters <- pl
    theta <- parhat[!obj$env$lrandom()]  ## ALL top-level parameters
    b <- parhat[obj$env$lrandom()]
    if (!is.null(obj$env$spHess)) {
        Hbb <- obj$env$spHess(parhat, random=TRUE) ## Needed later for RE models
    }
    ## #################################################################
    ## 1. Get partial derivatives of theta_hat wrt to yobs
    ## Note: length(yobs) much greater that length(theta)
    ##       ==> Reverse mode AD is suitable !
    ## #################################################################
    ## Move 'yobs' from data -> parameters (preserve C++ template order!)
    obj$env$parameters <- c(list(yobs = yobs), obj$env$parameters)
    obj$env$data$yobs <- NULL
    obj$retape()
    ## New objective parameters: (yobs, b, theta)
    nobs <- length(obj$env$parameters$yobs)
    nb <- length(obj$env$random)
    ntheta <- length(obj$env$par) - nobs - nb
    TMB::config(tmbad.atomic_sparse_log_determinant=0, DLL="RTMB") ## TMB FIXME
    F <- RTMB::GetTape(obj)
    r <- obj$env$random ## Which are random
    p <- tail(1:(nobs+ntheta), ntheta) ## Which are parameters *after* removing random
    ThetaHat <- F$laplace(r)$newton(p)
    J <- ThetaHat$jacobian(ThetaHat$par())
    ## Extra stuff we need in (3)
    F. <- F$jacfun() ## (yobs, [b], theta) -> (yobs, [b], theta)
    F. <- RTMB::MakeTape(function(y) F.( c(y, parhat) ) [r] , yobs) ## yobs -> b
    Hby <- F.$jacfun(sparse=TRUE)(yobs)
    ## #################################################################
    ## 2. Get partial derivatives of mu_hat wrt to theta for *fixed* yobs
    ## Note: length(mu) much greater that length(theta)
    ##       ==> Forward mode AD is suitable !
    ## #################################################################
    obj$env$data$yobs <- yobs
    obj$env$parameters$yobs <- NULL
    obj$retape()
    r <- obj$env$random ## Which are now random
    F <- RTMB::GetTape(obj)
    Bhat <- F$newton(r) ## theta -> bhat
    obj$env$data$doPredict <- as.double(1) ## Enable prediction of 'mu'
    obj$env$data$whichPredict <- as.double(1:nobs)
    obj$env$ADreport <- TRUE ## Modify return value from Cpp
    obj$retape()
    F <- RTMB::GetTape(obj) ## (b, theta) -> mu
    ## This doesn't work:
    ## MuHat <- MakeTape(function(theta)F(c(Bhat(theta), theta)), theta)
    MuHat <- RTMB::MakeTape(function(theta) {
        par <- RTMB::advector(nb + ntheta) ## glmmTMB mixes order of parameters and random effects...
        r <- obj$env$lrandom()
        par[r] <- Bhat(theta)
        par[!r] <- theta
        F(par)
    } , theta)
    ## 'Adjoint trick'
    T2 <- RTMB::MakeTape(function(weight) {
        WT <- RTMB::MakeTape(function(th) sum(MuHat(th) * weight), theta)
        WT$jacfun()(RTMB::advector(theta))
    }, rep(1,nobs))
    J2 <- T2$jacobian(yobs)
    if (diag) {
        term1 <- colSums(J*J2)
    } else {
        term1 <- t(J) %*% J2
    }
    ## #################################################################
    ## 3. Get partial derivatives of mu_hat wrt yobs for fixed theta
    ## Note: Tricky!
    ##       
    ## #################################################################
    term2 <- 0
    if (has.random) {
        F2 <- RTMB::MakeTape(function(b) {
            par <- RTMB::advector(nb + ntheta) ## glmmTMB mixes order of parameters and random effects...
            r <- obj$env$lrandom()
            par[r] <- b
            par[!r] <- theta
            F(par)
        }, b) ## (b) -> mu
        Hmb <- F2$jacfun(sparse=TRUE)(b) ## sparse deriv mu wrt b
        ## Implicit function theorem gives final partial deriv matrix:
        ##   - Hby %*% solve(Hbb) %*% Hbm
        ## of which we need the diagonal.
        ## Because mu and yobs link to the same random effects, all required b-cliques are part of Hbb !
        ## It follows that we can replace solve(Hbb) by its subset iH !
        if (diag) {
            iH <- TMB:::solveSubset(Hbb)
            ## FIXME: not quite sure why this coercion is needed?
            ## Error: no more error handlers available (recursive errors?); invoking 'abort' restart
            tmpM <- Hby * (  iH %*% t(Hmb) )
            term2 <- -colSums( as(tmpM, "matrix") )
        } else {
            term2 <- -t(Hby) %*% solve(Hbb, t(Hmb))
        }
    }
    term1 + term2
}

#' @export
rstudent.glmmTMB <- function (model, ...) {
    ## FIXME: more careful about 0-weight, na.exclude, naming, etc.
    check_dots(...)
    rp <- residuals(model, type = "pearson")
    r <- residuals(model, type = "deviance")
    hatval <- hatvalues(model)
    r <- sign(r) * sqrt(r^2 + (hatvalues * rp^2)/(1 - hatvalues))
    r[is.infinite(r)] <- NaN
    ## FIXME: account for infl$sigma
    ## (sd when obs i is dropped); should divide
    r
}

